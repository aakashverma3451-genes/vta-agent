"""Phase 9D — parent prodrug vs active triphosphate docking at the HCV NS5B catalytic site.

Phase 1 left an open question: when nucleotide-analog antivirals are docked as their PARENT
prodrugs (because curated triphosphate forms were not committed), how much of a poor recovery
is the prodrug artifact rather than a genuine scoring limitation? This script answers it
without fabricating anything: the HCV NS5B NI active-site set already contains BOTH real
parent/prodrug forms AND ChEMBL-curated triphosphate active forms, so we dock both groups
into the real NS5B catalytic site (2XI3 chain A, with the experimental catalytic Mg2+
retained) and compare Vina affinity (ΔG) distributions side by side.

Honesty:
  * Metal is read from the experimental 2XI3 structure (curated, not fabricated).
  * Ligands are docked AS GIVEN — parents as parents, triphosphates as triphosphates — with
    no species resolution, because the parent-vs-active-form difference is the measurement.
  * TiLV PB1 active-form re-dock is a labelled skip: its controls have no curated
    triphosphate SMILES (vta/chem/species.py active_form_smiles=None) and we will not
    fabricate them.
"""
from __future__ import annotations

import argparse
import json
import os
import statistics
import subprocess
import sys
from pathlib import Path

import requests

from vta.nodes.docking import _prep_ligand, _run_vina, _vina_bin

ACTIVES = Path("vta/data/phase9_actives/hcv_ns5b_ni_active_site.json")
OUT = Path("outputs/phase9/phase9d_active_form_comparison.json")
PDB_ID = "2XI3"             # HCV NS5B genotype 1a + GTP + catalytic Mg2+
CHAIN = "A"
CATALYTIC_ASP = ("220", "318", "319")   # GDD motif
CATALYTIC_MG_RESNUMS = {"1004", "1005"}  # the two active-site Mg (1006 is distal/surface)
TRIPHOSPHATE_MOTIF = "P(=O)(O)OP(=O)(O)OP(=O)(O)O"
_LIG_CACHE = "structures/.ligand_pdbqt_9d"


def fetch_pdb(pdb_id: str) -> str:
    r = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=90)
    r.raise_for_status()
    return r.text


def gdd_center(pdb_text: str, chain: str = CHAIN) -> list[float]:
    """Box center = centroid of the catalytic GDD aspartate CA atoms (frame-correct)."""
    pts = []
    for l in pdb_text.splitlines():
        if (l.startswith("ATOM") and l[21] == chain and l[12:16].strip() == "CA"
                and l[22:26].strip() in CATALYTIC_ASP and l[17:20].strip() == "ASP"):
            pts.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    n = len(pts)
    return [round(sum(p[i] for p in pts) / n, 3) for i in range(3)]


def build_receptor(pdb_text: str, out_pdbqt: str) -> str | None:
    """Chain-A protein + experimental catalytic Mg2+ → rigid receptor PDBQT (OpenBabel)."""
    keep = []
    for l in pdb_text.splitlines():
        if l.startswith("ATOM") and l[21] == CHAIN:
            keep.append(l)
        elif (l.startswith("HETATM") and l[21] == CHAIN and l[17:20].strip() == "MG"
              and l[22:26].strip() in CATALYTIC_MG_RESNUMS):
            keep.append(l)   # retain the curated catalytic metal
    if not keep:
        return None
    rec_pdb = out_pdbqt + ".receptor.pdb"
    Path(rec_pdb).write_text("\n".join(keep) + "\nEND\n")
    try:
        from openbabel import pybel
        mol = next(pybel.readfile("pdb", rec_pdb))
        mol.addh()
        mol.write("pdbqt", out_pdbqt, opt={"r": True}, overwrite=True)
    except Exception:
        return None
    return out_pdbqt if os.path.exists(out_pdbqt) and os.path.getsize(out_pdbqt) > 0 else None


def _largest_fragment(smiles: str) -> str:
    if "." not in smiles:
        return smiles
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
        return Chem.MolToSmiles(max(frags, key=lambda m: m.GetNumHeavyAtoms()))
    except Exception:
        return max(smiles.split("."), key=len)


def is_triphosphate(smiles: str) -> bool:
    return TRIPHOSPHATE_MOTIF in (smiles or "")


def run(limit: int | None = None) -> dict:
    vina = _vina_bin()
    if not vina:
        raise SystemExit("Vina not available")
    actives = json.loads(ACTIVES.read_text())["actives"]
    if limit:
        # keep a mix of both groups for a smoke run
        tri = [a for a in actives if is_triphosphate(a["canonical_smiles"])][:limit]
        par = [a for a in actives if not is_triphosphate(a["canonical_smiles"])][:limit]
        actives = tri + par

    pdb_text = fetch_pdb(PDB_ID)
    center = gdd_center(pdb_text)
    os.makedirs("structures", exist_ok=True)
    os.makedirs(_LIG_CACHE, exist_ok=True)
    receptor = build_receptor(pdb_text, f"structures/ns5b_{PDB_ID}_rec.pdbqt")
    if not receptor:
        raise SystemExit("NS5B receptor prep failed")

    results = []
    for a in actives:
        smi = _largest_fragment(a["canonical_smiles"])
        cid = a["molecule_chembl_id"]
        lig = {"smiles": smi, "chembl_id": cid, "dock_cache_id": f"9d_{cid}"}
        try:
            lig_pdbqt = _prep_ligand(lig, _LIG_CACHE)
        except Exception:
            lig_pdbqt = None
        dG = None
        if lig_pdbqt:
            out = f"structures/9d_{cid}_pose.pdbqt"
            try:
                dG = _run_vina(vina, receptor, lig_pdbqt, center, out)
            except Exception:
                dG = None
        results.append({
            "molecule_chembl_id": cid,
            "pref_name": a.get("pref_name"),
            "form": "triphosphate" if is_triphosphate(a["canonical_smiles"]) else "parent_or_prodrug",
            "dG": dG,
            "docked": dG is not None,
        })

    def _stats(group: str) -> dict:
        dgs = [r["dG"] for r in results if r["form"] == group and r["dG"] is not None]
        if not dgs:
            return {"n_docked": 0, "median_dG": None, "best_dG": None, "worst_dG": None}
        return {"n_docked": len(dgs), "median_dG": round(statistics.median(dgs), 3),
                "best_dG": round(min(dgs), 3), "worst_dG": round(max(dgs), 3)}

    tri_stats = _stats("triphosphate")
    par_stats = _stats("parent_or_prodrug")

    # Matched pair: sofosbuvir (parent) vs its triphosphate GS-461203, if both docked.
    by_name = {r["pref_name"]: r for r in results}
    matched = None
    if by_name.get("SOFOSBUVIR") and by_name.get("GS-461203"):
        matched = {
            "parent": {"name": "SOFOSBUVIR", "dG": by_name["SOFOSBUVIR"]["dG"]},
            "active_form": {"name": "GS-461203 (sofosbuvir triphosphate)",
                            "dG": by_name["GS-461203"]["dG"]},
        }

    delta = (None if tri_stats["median_dG"] is None or par_stats["median_dG"] is None
             else round(tri_stats["median_dG"] - par_stats["median_dG"], 3))
    payload = {
        "phase": "9D parent-vs-active-form",
        "target": "HCV NS5B catalytic site",
        "structure": f"{PDB_ID} chain {CHAIN} (GTP + catalytic Mg2+)",
        "metal": {"present": True, "ion": "Mg2+", "resnums": sorted(CATALYTIC_MG_RESNUMS),
                  "provenance": f"experimental {PDB_ID} HETATM (curated, not fabricated)"},
        "box_center": center,
        "box_center_source": "centroid of catalytic GDD aspartates 220/318/319 CA",
        "scoring": "AutoDock Vina ΔG (kcal/mol); ligands docked AS GIVEN (no species resolution)",
        "n_actives": len(actives),
        "triphosphate_active_form": tri_stats,
        "parent_or_prodrug": par_stats,
        "median_dG_triphosphate_minus_parent": delta,
        "matched_pair_sofosbuvir": matched,
        "tilv_active_form": {
            "status": "labelled skip",
            "reason": "TiLV PB1 controls have no curated triphosphate active-form SMILES "
                      "(vta/chem/species.py active_form_smiles=None); not fabricated.",
        },
        "interpretation": (
            "More negative ΔG = stronger predicted binding. If the triphosphate group docks "
            "markedly stronger than the parent/prodrug group at the metal-containing catalytic "
            "site, the poor parent recovery is substantially a prodrug-input artifact, not "
            "purely a scoring limitation. ΔG is not enrichment: HCV NS5B NI is un-benchmarkable "
            "by matched decoys (Phase 9B-3 meta-finding), so this is a within-set score "
            "comparison, reported as such."
        ),
        "results": results,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2))
    _write_md(payload)
    return payload


def _write_md(p: dict) -> None:
    tri, par = p["triphosphate_active_form"], p["parent_or_prodrug"]
    md = [
        "# Phase 9D: Parent prodrug vs active triphosphate (HCV NS5B catalytic site)",
        "",
        f"Structure: {p['structure']}. Box center {p['box_center']} "
        f"({p['box_center_source']}). Metal: {p['metal']['ion']} retained from "
        f"{p['metal']['provenance']}.",
        f"Scoring: {p['scoring']}.",
        "",
        "## ΔG by form (more negative = stronger predicted binding)",
        "",
        "| Form | n docked | median ΔG | best ΔG | worst ΔG |",
        "|---|---|---|---|---|",
        f"| Triphosphate active form | {tri['n_docked']} | {tri['median_dG']} | {tri['best_dG']} | {tri['worst_dG']} |",
        f"| Parent / prodrug | {par['n_docked']} | {par['median_dG']} | {par['best_dG']} | {par['worst_dG']} |",
        "",
        f"Median ΔG (triphosphate − parent): **{p['median_dG_triphosphate_minus_parent']}** "
        "kcal/mol (negative = triphosphate binds stronger).",
        "",
    ]
    if p.get("matched_pair_sofosbuvir"):
        m = p["matched_pair_sofosbuvir"]
        md += [
            "## Matched pair",
            "",
            f"- Sofosbuvir (parent): ΔG {m['parent']['dG']}",
            f"- {m['active_form']['name']}: ΔG {m['active_form']['dG']}",
            "",
        ]
    md += [
        "## Interpretation",
        "",
        p["interpretation"],
        "",
        "## TiLV active-form",
        "",
        f"- {p['tilv_active_form']['status']}: {p['tilv_active_form']['reason']}",
    ]
    OUT.with_suffix(".md").write_text("\n".join(md) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None, help="smoke run: N per form")
    args = ap.parse_args()
    p = run(limit=args.limit)
    print(f"triphosphate: {p['triphosphate_active_form']}")
    print(f"parent/prodrug: {p['parent_or_prodrug']}")
    print(f"median ΔG (tri - parent) = {p['median_dG_triphosphate_minus_parent']}")
    print(OUT)


if __name__ == "__main__":
    main()
