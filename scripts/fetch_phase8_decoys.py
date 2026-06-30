"""Fetch Phase 8 property-matched presumed decoys from ChEMBL."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import requests

from vta.data.ligands import load_ligands
from vta.eval.splits import scaffold_key

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
OUT_SMI = Path("vta/data/decoys_cache/phase8_tilv_pb1_chembl_matched.smi")
OUT_JSON = Path("outputs/phase8/phase8_decoy_provenance.json")


def _props(smiles: str) -> dict:
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    return {
        "mw": Descriptors.MolWt(mol),
        "alogp": Crippen.MolLogP(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "hbd": Lipinski.NumHDonors(mol),
        "rtb": Lipinski.NumRotatableBonds(mol),
        "scaffold": scaffold_key(smiles),
    }


def _query(active: dict, limit: int, pages: int) -> list[dict]:
    p = _props(active["smiles"])
    base = {
        "limit": min(80, max(20, limit * 2)),
        "molecule_properties__mw_freebase__gte": max(100, int(p["mw"] - 40)),
        "molecule_properties__mw_freebase__lte": int(p["mw"] + 40),
        "format": "json",
    }
    response = None
    rows = []
    for width in (40, 80, 120):
        try:
            for page in range(pages):
                params = {
                    **base,
                    "offset": page * base["limit"],
                    "molecule_properties__mw_freebase__gte": max(100, int(p["mw"] - width)),
                    "molecule_properties__mw_freebase__lte": int(p["mw"] + width),
                }
                response = requests.get(CHEMBL, params=params, timeout=30)
                response.raise_for_status()
                for mol in response.json().get("molecules", []):
                    structures = mol.get("molecule_structures") or {}
                    smiles = structures.get("canonical_smiles")
                    chembl_id = mol.get("molecule_chembl_id")
                    if smiles and chembl_id and _locally_matched(p, smiles):
                        rows.append({
                            "name": mol.get("pref_name") or chembl_id,
                            "chembl_id": chembl_id,
                            "smiles": smiles,
                            "matched_active": active["name"],
                        })
                if len(rows) >= limit:
                    return rows
            if rows:
                return rows
        except requests.RequestException:
            if width == 120:
                raise
            time.sleep(1)
    return rows


def _locally_matched(active_props: dict, smiles: str) -> bool:
    props = _props(smiles)
    if not props:
        return False
    return (
        abs(props["mw"] - active_props["mw"]) <= 120
        and abs(props["alogp"] - active_props["alogp"]) <= 2.0
        and abs(props["hba"] - active_props["hba"]) <= 4
        and abs(props["hbd"] - active_props["hbd"]) <= 4
        and abs(props["rtb"] - active_props["rtb"]) <= 6
    )


def fetch(per_active: int = 30, delay: float = 0.5, pages: int = 1) -> dict:
    actives = [lig for lig in load_ligands() if lig.get("positive_control")]
    active_ids = {lig.get("chembl_id") for lig in load_ligands()}
    active_scaffolds = {scaffold_key(lig["smiles"]) for lig in actives}
    selected: list[dict] = []
    seen = set(active_ids)

    for active in actives:
        kept = 0
        try:
            candidates = _query(active, per_active, pages)
        except requests.RequestException as exc:
            selected.append({
                "name": f"FETCH_FAILED_{active['name']}",
                "chembl_id": None,
                "smiles": None,
                "matched_active": active["name"],
                "error": f"{type(exc).__name__}: {exc}",
            })
            continue
        for row in candidates:
            if row["chembl_id"] in seen:
                continue
            if scaffold_key(row["smiles"]) in active_scaffolds:
                continue
            row["property_match_method"] = "ChEMBL MW/logP/HBA/HBD window; scaffold-distinct"
            selected.append(row)
            seen.add(row["chembl_id"])
            kept += 1
            if kept >= per_active:
                break
        time.sleep(delay)

    OUT_SMI.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    usable = [row for row in selected if row.get("smiles") and row.get("chembl_id")]
    existing_n = sum(1 for line in OUT_SMI.read_text().splitlines() if line.strip()) if OUT_SMI.exists() else 0
    preserved_existing = existing_n > len(usable)
    if not preserved_existing:
        OUT_SMI.write_text("\n".join(
            f"{row['smiles']}\t{row['name']}\t{row['chembl_id']}\tmatched_to={row['matched_active']}"
            for row in usable
        ) + "\n")
    provenance = {
        "source": "ChEMBL molecule API",
        "method": "property-matched presumed decoys; MW/logP/HBA/HBD windows; scaffold-distinct",
        "target": "TiLV PB1",
        "per_active_requested": per_active,
        "pages_requested": pages,
        "actives": [a["name"] for a in actives],
        "n_decoys": len(usable),
        "existing_decoys_preserved": preserved_existing,
        "n_failures": len(selected) - len(usable),
        "output_smi": str(OUT_SMI),
        "caveat": "Presumed decoys are not experimentally verified inactives.",
        "decoys": usable,
        "failures": [row for row in selected if not row.get("smiles")],
    }
    OUT_JSON.write_text(json.dumps(provenance, indent=2))
    return provenance


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--per-active", type=int, default=30)
    parser.add_argument("--pages", type=int, default=1)
    args = parser.parse_args()
    out = fetch(args.per_active, pages=args.pages)
    print(f"fetched {out['n_decoys']} ChEMBL property-matched presumed decoys")
    print(out["output_smi"])


if __name__ == "__main__":
    main()
