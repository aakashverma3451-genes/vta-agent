"""Fetch versioned Phase 9 known-actives from ChEMBL assays."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import requests

CHEMBL_TARGET = "https://www.ebi.ac.uk/chembl/api/data/target.json"
CHEMBL_ACTIVITY = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
OUT = Path("vta/data/phase9_actives")


def search_targets(query: str, limit: int = 10) -> list[dict]:
    r = requests.get(
        CHEMBL_TARGET,
        params={"pref_name__icontains": query, "limit": limit, "format": "json"},
        timeout=30,
    )
    r.raise_for_status()
    return r.json().get("targets", [])


def fetch_activities(target_chembl_id: str, limit: int = 1000, pages: int = 5) -> list[dict]:
    rows = []
    for page in range(pages):
        params = {
            "target_chembl_id": target_chembl_id,
            "standard_type__in": "IC50,EC50,Ki,Kd",
            "standard_units": "nM",
            "standard_value__isnull": "false",
            "pchembl_value__isnull": "false",
            "limit": limit,
            "offset": page * limit,
            "format": "json",
        }
        r = requests.get(CHEMBL_ACTIVITY, params=params, timeout=45)
        r.raise_for_status()
        page_rows = r.json().get("activities", [])
        rows.extend(page_rows)
        if len(page_rows) < limit:
            break
        time.sleep(0.5)
    return rows


def collapse_actives(rows: list[dict], max_nm: float = 10000.0) -> list[dict]:
    by_mol: dict[str, dict] = {}
    for row in rows:
        mol_id = row.get("molecule_chembl_id")
        try:
            value = float(row.get("standard_value"))
            pchembl = float(row.get("pchembl_value"))
        except (TypeError, ValueError):
            continue
        relation = row.get("standard_relation") or "="
        if not mol_id or value > max_nm or relation.strip() not in {"=", "<", "<="}:
            continue
        current = by_mol.get(mol_id)
        if current is None or pchembl > current["best_pchembl"]:
            by_mol[mol_id] = {
                "molecule_chembl_id": mol_id,
                "canonical_smiles": row.get("canonical_smiles"),
                "pref_name": row.get("molecule_pref_name") or mol_id,
                "best_pchembl": pchembl,
                "activity_nm": value,
                "activity_type": row.get("standard_type"),
                "assay_chembl_id": row.get("assay_chembl_id"),
                "document_chembl_id": row.get("document_chembl_id"),
                "target_chembl_id": row.get("target_chembl_id"),
                "provenance": "ChEMBL activity endpoint",
            }
    return sorted(by_mol.values(), key=lambda r: r["best_pchembl"], reverse=True)


def write_actives(label: str, target: dict, actives: list[dict]) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = {
        "label": label,
        "target": target,
        "n_actives": len(actives),
        "activity_cutoff_nm": 10000,
        "actives": actives,
        "caveat": "ChEMBL assay actives require assay-type review before publication claims.",
    }
    path = OUT / f"{label}.json"
    path.write_text(json.dumps(payload, indent=2))
    return path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--query", default="HCV NS5B polymerase")
    parser.add_argument("--target-id")
    parser.add_argument("--label", default="hcv_ns5b")
    parser.add_argument("--pages", type=int, default=5)
    args = parser.parse_args()

    if args.target_id:
        target = {"target_chembl_id": args.target_id, "pref_name": args.query}
    else:
        targets = search_targets(args.query)
        if not targets:
            raise SystemExit(f"no ChEMBL targets found for {args.query!r}")
        target = targets[0]
    rows = fetch_activities(target["target_chembl_id"], pages=args.pages)
    actives = collapse_actives(rows)
    path = write_actives(args.label, target, actives)
    print(f"target={target.get('target_chembl_id')} {target.get('pref_name')}")
    print(f"actives={len(actives)}")
    print(path)


if __name__ == "__main__":
    main()
