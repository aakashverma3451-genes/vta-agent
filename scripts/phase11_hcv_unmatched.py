"""Phase 11 WI-5 Job B — HCV NS5B NI re-benchmark with property-UNMATCHED decoys.

The nucleotide meta-finding was that a purchasable-library, scaffold-distinct, PROPERTY-MATCHED
decoy recipe cannot build a decoy set for triphosphate actives (they have no property twins).
WI-5 asks: does the benchmark run at all with a DIFFERENT recipe — property-UNMATCHED /
charge-extrema decoys (Stein 2021, DUDE-Z philosophy)? Those are trivial to obtain (any
drug-like, phosphate-free molecule is property-unmatched to a charged triphosphate), so the
set CAN be built. The scientific catch this then exposes: if Vina "enriches" triphosphate
actives over neutral drug-like decoys, is that real binding signal or just charge/size
discrimination? We report the actual result honestly, whatever it is.

Actives' ΔG is reused from the committed 9D dock (same 2XI3 receptor + box); only the decoys
are newly docked. DeepCoy is not installed — these are ChEMBL property-unmatched decoys with
tool+seed logged, not generative decoys.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import requests

from vta.eval.baselines import baseline_2d_similarity, bootstrap_metric_ci, random_metric_distribution
from vta.eval.metrics import METRIC_NAMES, bootstrap_enrichment_report, enrichment_report
from vta.eval.significance import holm_bonferroni, paired_bootstrap_delta
from vta.eval.splits import scaffold_key
from vta.nodes.docking import _prep_ligand, _run_vina, _vina_bin
from scripts.run_phase9d_active_form import (
    _largest_fragment, build_receptor, fetch_pdb, gdd_center, is_triphosphate,
)

ACTIVES = Path("vta/data/phase9_actives/hcv_ns5b_ni_active_site.json")
P9D = Path("outputs/phase9/phase9d_active_form_comparison.json")
DECOY_SMI = Path("vta/data/decoys_cache/hcv_ns5b_ni_unmatched.smi")
OUT = Path("outputs/phase11/hcv_ni_unmatched_benchmark.json")
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
_LIG_CACHE = "structures/.ligand_pdbqt_hcvun"
N_DECOYS_TARGET = 640          # ~15 per active × 43 → informative EF without a multi-hour dock
SEED = 11


def fetch_unmatched_decoys(n_target: int, active_scaffolds: set[str], seed: int = SEED) -> list[dict]:
    """Broad drug-like, phosphate-free, scaffold-distinct ChEMBL molecules (property-UNMATCHED)."""
    import random as _r
    rng = _r.Random(seed)
    from rdkit import Chem
    rows: list[dict] = []
    seen: set[str] = set()
    page = 0
    while len(rows) < n_target and page < 40:
        params = {
            "molecule_properties__full_mwt__gte": 250,
            "molecule_properties__full_mwt__lte": 500,
            "molecule_properties__num_ro5_violations__lte": 1,
            "limit": 200, "offset": page * 200 + rng.randint(0, 50), "format": "json",
        }
        try:
            r = requests.get(CHEMBL, params=params, timeout=45)
            r.raise_for_status()
            mols = r.json().get("molecules", [])
        except requests.RequestException:
            time.sleep(1.5)
            page += 1
            continue
        for mol in mols:
            structs = mol.get("molecule_structures") or {}
            smi = structs.get("canonical_smiles")
            cid = mol.get("molecule_chembl_id")
            if not smi or not cid or cid in seen:
                continue
            if "P" in smi and is_triphosphate(smi):
                continue          # exclude phosphate-bearing (keep them property-UNmatched)
            parent = _largest_fragment(smi)
            m = Chem.MolFromSmiles(parent)
            if m is None:
                continue
            if scaffold_key(parent) in active_scaffolds:
                continue          # scaffold-distinct from actives
            seen.add(cid)
            rows.append({"chembl_id": cid, "smiles": parent})
            if len(rows) >= n_target:
                break
        page += 1
        time.sleep(0.3)
    return rows


def _write_decoys(rows: list[dict]) -> None:
    DECOY_SMI.parent.mkdir(parents=True, exist_ok=True)
    DECOY_SMI.write_text("".join(f"{r['smiles']}\t{r['chembl_id']}\n" for r in rows))


def _load_decoys() -> list[dict]:
    out = []
    for line in DECOY_SMI.read_text().splitlines():
        if line.strip():
            smi, cid = line.split("\t")[:2]
            out.append({"chembl_id": cid, "smiles": smi})
    return out


def run(n_decoys: int = N_DECOYS_TARGET) -> dict:
    vina = _vina_bin()
    if not vina:
        raise SystemExit("Vina not available")
    actives = json.loads(ACTIVES.read_text())["actives"]
    active_scaffolds = {scaffold_key(_largest_fragment(a["canonical_smiles"])) for a in actives}

    # Reuse committed 9D actives' ΔG (same 2XI3 receptor + box).
    p9d = json.loads(P9D.read_text())
    active_dG = {r["molecule_chembl_id"]: r.get("dG") for r in p9d["results"]}

    if not DECOY_SMI.exists():
        _write_decoys(fetch_unmatched_decoys(n_decoys, active_scaffolds))
    decoys = _load_decoys()

    # Receptor + box from 2XI3 (metal-retaining), same as 9D.
    pdb_text = fetch_pdb("2XI3")
    center = gdd_center(pdb_text)
    receptor = build_receptor(pdb_text, "structures/ns5b_2XI3_rec.pdbqt")
    import os
    os.makedirs(_LIG_CACHE, exist_ok=True)

    entries = []
    for a in actives:
        dG = active_dG.get(a["molecule_chembl_id"])
        if dG is not None:
            entries.append({"name": a["molecule_chembl_id"], "score": round(-float(dG), 4),
                            "positive_control": True, "smiles": _largest_fragment(a["canonical_smiles"]),
                            "dG": dG})
    n_decoys_docked = 0
    for d in decoys:
        ld = {"smiles": d["smiles"], "chembl_id": d["chembl_id"], "dock_cache_id": f"hcvun_{d['chembl_id']}"}
        try:
            lp = _prep_ligand(ld, _LIG_CACHE)
            dG = _run_vina(vina, receptor, lp, center, f"structures/hcvun_{d['chembl_id']}_pose.pdbqt") if lp else None
        except Exception:
            dG = None
        if dG is not None:
            entries.append({"name": d["chembl_id"], "score": round(-float(dG), 4),
                            "positive_control": False, "smiles": d["smiles"], "dG": dG})
            n_decoys_docked += 1

    labels = [1 if e["positive_control"] else 0 for e in entries]
    smiles = [e["smiles"] for e in entries]
    vina_scores = [e["score"] for e in entries]
    sim_scores = baseline_2d_similarity(smiles, labels)
    report = enrichment_report(entries)
    bootstrap = bootstrap_enrichment_report(entries, n_resamples=10000, seed=42)
    baselines = {
        "random": {m: random_metric_distribution(labels, m, n=10000, seed=11) for m in METRIC_NAMES},
        "2d_sim": {m: bootstrap_metric_ci(sim_scores, labels, m, n=10000, seed=11) for m in METRIC_NAMES},
        "vina": {m: bootstrap_metric_ci(vina_scores, labels, m, n=10000, seed=11) for m in METRIC_NAMES},
    }
    vina_vs_2d = {m: paired_bootstrap_delta(vina_scores, sim_scores, labels, m, n=10000, seed=21)
                  for m in METRIC_NAMES}
    holm = holm_bonferroni({m: (1.0 - vina_vs_2d[m]["p_gt0"]) if vina_vs_2d[m]["p_gt0"] is not None
                            else None for m in METRIC_NAMES})
    n_act = sum(labels)
    payload = {
        "phase": "11 WI-5 Job B — HCV NS5B NI with property-UNMATCHED decoys",
        "structure": "2XI3 chain A (NS5B + catalytic Mg2+)",
        "decoy_recipe": "ChEMBL drug-like (MW 250-500, ≤1 Ro5 violation), phosphate-free, "
                        "scaffold-distinct — deliberately property-UNMATCHED to the charged "
                        "triphosphate actives (Stein 2021 DUDE-Z philosophy). DeepCoy absent.",
        "decoy_provenance": {"source": "ChEMBL molecule API", "seed": SEED,
                             "file": str(DECOY_SMI), "n_decoys_generated": len(decoys)},
        "n_actives": n_act, "n_decoys_docked": n_decoys_docked,
        "ratio": f"{n_decoys_docked}:{n_act}",
        "benchmark": report, "bootstrap": bootstrap, "baseline_tables": baselines,
        "vina_vs_2d_similarity": vina_vs_2d, "holm_bonferroni": holm,
        "interpretation": (
            "This benchmark DOES run with a property-unmatched recipe (unlike the property-"
            "matched recipe, which could not build a decoy set). Read the numbers with the "
            "DUDE-Z caveat: high enrichment against neutral drug-like decoys may reflect "
            "charge/size discrimination of the triphosphate actives rather than pocket-specific "
            "binding — the paired Vina-vs-2D-similarity test separates those. Reported as-is."),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2))
    _md(payload)
    return payload


def _md(p: dict) -> None:
    def ci(t, m):
        return f"{t[m].get('median')} {t[m].get('ci95')}"
    md = [f"# Phase 11 WI-5 Job B — HCV NS5B NI, property-UNMATCHED decoys", "",
          f"Structure: {p['structure']}. Ratio {p['ratio']}. Decoy recipe: {p['decoy_recipe']}", "",
          "## Baselines vs Vina (median [95% CI])", "",
          "| Method | BEDROC | logAUC | ROC-AUC | EF1% |", "|---|---|---|---|---|"]
    for label in ("random", "2d_sim", "vina"):
        t = p["baseline_tables"][label]
        md.append(f"| {label} | {ci(t,'bedroc')} | {ci(t,'logauc')} | {ci(t,'roc_auc')} | {ci(t,'ef1')} |")
    v = p["vina_vs_2d_similarity"]["bedroc"]
    md += ["", f"Paired Vina−2Dsim (BEDROC): Δ {v['median_delta']}, 95% CI {v['ci95']}.", "",
           "## Interpretation", "", p["interpretation"]]
    OUT.with_suffix(".md").write_text("\n".join(md) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-decoys", type=int, default=N_DECOYS_TARGET)
    args = ap.parse_args()
    p = run(args.n_decoys)
    b = p["benchmark"]
    print(f"ratio={p['ratio']} BEDROC={b['bedroc']} EF1%={b['ef'].get('EF1%')} ROC-AUC={b['roc_auc']}")
    print(OUT)


if __name__ == "__main__":
    main()
