"""Assemble a SARS-CoV-2 Mpro dataset of experimentally-measured actives AND inactives.

Track B of Phase 9C needs a drug-like target with *real measured inactives* — the gold
standard that removes the presumed-decoy weakness entirely. Two sources, both behind
injectable seams so tests stay hermetic:

  * PRIMARY — COVID Moonshot (``postera-ai/COVID_moonshot_submissions``,
    ``covid_submissions_all_info.csv``): fluorescence IC50 (``f_avg_IC50`` µM) plus
    %-inhibition-at-50µM and explicit covalent annotations. Compounds assayed and found
    *inactive* are real measured inactives — exactly what a powered enrichment control set
    should be.
  * TOP-UP — ChEMBL target ``CHEMBL4523582`` (SARS-CoV-2 Mpro), to cross-check and expand
    the active count.

Activity labelling (not binding mode — that is stratified separately, see
``stratify_mpro_binding_mode.py``):
  * active   — measured IC50 ≤ ACTIVE_MAX_UM.
  * inactive — measured IC50 ≥ INACTIVE_MIN_UM, or assayed with %inhibition@50µM below
               INACTIVE_INHIB_PCT (assayed, did not inhibit).
  * else     — intermediate; kept out of the benchmark and counted, never silently pooled.

No SMILES or activity value is fabricated; rows without a measured value are dropped.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import time
from pathlib import Path
from typing import Callable, Optional

import requests

MOONSHOT_CSV = (
    "https://raw.githubusercontent.com/postera-ai/COVID_moonshot_submissions/"
    "master/covid_submissions_all_info.csv"
)
CHEMBL_ACTIVITY = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
CHEMBL_MPRO_TARGET = "CHEMBL4523582"  # SARS-CoV-2 "Replicase polyprotein 1ab" (Mpro)

OUT = Path("vta/data/mpro/mpro_dataset.json")

ACTIVE_MAX_UM = 10.0      # ≤ 10 µM measured IC50 → active
INACTIVE_MIN_UM = 50.0    # ≥ 50 µM measured IC50 → inactive
INACTIVE_INHIB_PCT = 25.0  # < 25% inhibition @ 50 µM → assayed-inactive
PUBLICATION_GRADE_PER_ACTIVE = 30

MoonshotFetch = Callable[[], str]
ChemblFetch = Callable[[int], list[dict]]


# ── network seams (monkeypatched in tests) ───────────────────────────────────
def default_moonshot_fetch(timeout: int = 120) -> str:
    r = requests.get(MOONSHOT_CSV, timeout=timeout)
    r.raise_for_status()
    return r.text


def default_chembl_fetch(pages: int = 6, limit: int = 1000, timeout: int = 60,
                         retries: int = 4) -> list[dict]:
    rows: list[dict] = []
    for page in range(pages):
        params = {
            "target_chembl_id": CHEMBL_MPRO_TARGET,
            "standard_type__in": "IC50",
            "standard_units": "nM",
            "limit": limit,
            "offset": page * limit,
            "format": "json",
        }
        last_exc: Optional[Exception] = None
        for attempt in range(retries):
            try:
                r = requests.get(CHEMBL_ACTIVITY, params=params, timeout=timeout)
                r.raise_for_status()
                page_rows = r.json().get("activities", [])
                break
            except requests.RequestException as exc:
                last_exc = exc
                time.sleep(1.5 * (attempt + 1))
        else:
            raise last_exc  # type: ignore[misc]
        rows.extend(page_rows)
        if len(page_rows) < limit:
            break
        time.sleep(0.4)
    return rows


# ── pure helpers ──────────────────────────────────────────────────────────────
def _num(value: object) -> Optional[float]:
    try:
        f = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return f if f == f else None  # drop NaN


def canonical_smiles(smiles: str | None) -> Optional[str]:
    """RDKit canonical SMILES (largest fragment) for dedup; falls back to raw string."""
    if not smiles:
        return None
    try:
        from rdkit import Chem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
        if frags:
            mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
        return Chem.MolToSmiles(mol)
    except Exception:
        return smiles


def parse_moonshot(csv_text: str,
                   active_max_um: float = ACTIVE_MAX_UM,
                   inactive_min_um: float = INACTIVE_MIN_UM,
                   inactive_inhib_pct: float = INACTIVE_INHIB_PCT) -> list[dict]:
    """Parse the Moonshot CSV into labelled, provenance-bearing records."""
    reader = csv.DictReader(io.StringIO(csv_text))
    records = []
    for row in reader:
        smiles = (row.get("SMILES") or "").strip()
        if not smiles:
            continue
        ic50 = _num(row.get("f_avg_IC50"))            # µM
        pic50 = _num(row.get("f_avg_pIC50"))
        inhib50 = _num(row.get("f_inhibition_at_50_uM"))
        if ic50 is not None and ic50 <= active_max_um:
            label = "active"
        elif (ic50 is not None and ic50 >= inactive_min_um) or (
                ic50 is None and inhib50 is not None and inhib50 < inactive_inhib_pct):
            label = "inactive"
        elif ic50 is not None or inhib50 is not None:
            label = "intermediate"
        else:
            continue  # no measurement at all → drop (not fabricated)
        covalent_annotation = (row.get("covalent_warhead")
                               or row.get("Covalent Fragment") or "").strip()
        records.append({
            "id": row.get("CID") or row.get("CDD_name"),
            "smiles": smiles,
            "label": label,
            "source": "moonshot",
            "measure": {"ic50_um": ic50, "pic50": pic50, "inhibition_at_50uM": inhib50},
            "covalent_annotation": covalent_annotation,
            "series": (row.get("series") or "").strip(),
        })
    return records


def parse_chembl(activities: list[dict],
                 active_pchembl: float = 6.0,
                 inactive_pchembl: float = 4.5) -> list[dict]:
    """Parse ChEMBL Mpro IC50 activities into labelled records (best per molecule)."""
    by_mol: dict[str, dict] = {}
    for a in activities:
        mol_id = a.get("molecule_chembl_id")
        smiles = a.get("canonical_smiles")
        if not mol_id or not smiles:
            continue
        value_nm = _num(a.get("standard_value"))
        pchembl = _num(a.get("pchembl_value"))
        comment = (a.get("activity_comment") or "").strip().lower()
        # Label from pChEMBL when present, else from an explicit inactive comment.
        if pchembl is not None and pchembl >= active_pchembl:
            label = "active"
        elif (pchembl is not None and pchembl <= inactive_pchembl) or comment in {
                "not active", "inactive", "no activity"}:
            label = "inactive"
        elif pchembl is not None:
            label = "intermediate"
        else:
            continue
        rank = pchembl if pchembl is not None else 0.0
        cur = by_mol.get(mol_id)
        if cur is None or rank > cur["_rank"]:
            by_mol[mol_id] = {
                "id": mol_id,
                "smiles": smiles,
                "label": label,
                "source": "chembl",
                "measure": {"ic50_nm": value_nm, "pchembl": pchembl,
                            "activity_comment": a.get("activity_comment")},
                "covalent_annotation": "",
                "series": "",
                "_rank": rank,
            }
    for rec in by_mol.values():
        rec.pop("_rank", None)
    return list(by_mol.values())


def assemble(moonshot_fetch: MoonshotFetch = default_moonshot_fetch,
             chembl_fetch: ChemblFetch = default_chembl_fetch,
             pages: int = 6) -> dict:
    """Merge Moonshot (primary) + ChEMBL (top-up), dedup by canonical SMILES."""
    provenance_sources = []
    moonshot_records: list[dict] = []
    try:
        moonshot_records = parse_moonshot(moonshot_fetch())
        provenance_sources.append({"source": "moonshot", "url": MOONSHOT_CSV,
                                   "n_parsed": len(moonshot_records)})
    except Exception as exc:  # degrade to ChEMBL-primary, labelled honestly
        provenance_sources.append({"source": "moonshot", "url": MOONSHOT_CSV,
                                   "error": f"{type(exc).__name__}: {exc}"})

    chembl_records: list[dict] = []
    try:
        chembl_records = parse_chembl(chembl_fetch(pages))
        provenance_sources.append({"source": "chembl", "target": CHEMBL_MPRO_TARGET,
                                   "n_parsed": len(chembl_records)})
    except Exception as exc:
        provenance_sources.append({"source": "chembl", "target": CHEMBL_MPRO_TARGET,
                                   "error": f"{type(exc).__name__}: {exc}"})

    # Moonshot first (preferred); ChEMBL only adds molecules Moonshot didn't cover.
    merged: list[dict] = []
    seen: set[str] = set()
    for rec in moonshot_records + chembl_records:
        key = canonical_smiles(rec["smiles"]) or rec["smiles"]
        if key in seen:
            continue
        seen.add(key)
        rec = {**rec, "canonical_smiles": key}
        merged.append(rec)

    actives = [r for r in merged if r["label"] == "active"]
    inactives = [r for r in merged if r["label"] == "inactive"]
    intermediate = [r for r in merged if r["label"] == "intermediate"]

    def _by_source(records: list[dict]) -> dict:
        out: dict[str, int] = {}
        for r in records:
            out[r["source"]] = out.get(r["source"], 0) + 1
        return out

    return {
        "target": "SARS-CoV-2 Mpro (3CLpro)",
        "sources": provenance_sources,
        "activity_thresholds": {
            "active_max_um": ACTIVE_MAX_UM,
            "inactive_min_um": INACTIVE_MIN_UM,
            "inactive_inhibition_at_50uM_pct": INACTIVE_INHIB_PCT,
            "chembl_active_pchembl": 6.0,
            "chembl_inactive_pchembl": 4.5,
        },
        "n_actives": len(actives),
        "n_inactives": len(inactives),
        "n_intermediate_excluded": len(intermediate),
        "actives_by_source": _by_source(actives),
        "inactives_by_source": _by_source(inactives),
        "caveat": (
            "Activity labels are experimentally measured (Moonshot fluorescence IC50 / "
            "%inhibition; ChEMBL pChEMBL). Binding mode (covalent vs non-covalent) is "
            "stratified separately before benchmarking; do not pool the two."
        ),
        "actives": actives,
        "inactives": inactives,
        "intermediate": intermediate,
    }


def write(dataset: dict) -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(dataset, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pages", type=int, default=6)
    args = parser.parse_args()
    dataset = assemble(pages=args.pages)
    write(dataset)
    print(f"actives={dataset['n_actives']} (by source {dataset['actives_by_source']})")
    print(f"inactives={dataset['n_inactives']} (by source {dataset['inactives_by_source']})")
    print(f"intermediate_excluded={dataset['n_intermediate_excluded']}")
    print(OUT)


if __name__ == "__main__":
    main()
