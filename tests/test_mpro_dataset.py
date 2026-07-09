"""Hermetic tests for the Mpro dataset assembly (Moonshot + ChEMBL seams faked)."""
from __future__ import annotations

from scripts.fetch_mpro_dataset import (
    assemble,
    canonical_smiles,
    parse_chembl,
    parse_moonshot,
)

# Minimal Moonshot-shaped CSV: an active, a high-IC50 inactive, a low-inhibition inactive,
# an intermediate, and a no-measurement row that must be dropped.
MOONSHOT_CSV = (
    "SMILES,CID,f_avg_IC50,f_avg_pIC50,f_inhibition_at_50_uM,covalent_warhead,series\n"
    "CCO,ACT-1,0.5,6.3,95,,aminopyridine\n"            # active (0.5 µM)
    "CCN,INA-1,80,4.1,,,quinolone\n"                    # inactive (80 µM)
    "CCC,INA-2,,,5,,benzo\n"                            # inactive (5% inhib @ 50µM)
    "CCCC,INT-1,25,4.6,40,,misc\n"                      # intermediate (25 µM)
    "CCCCC,DROP-1,,,,,\n"                               # no measurement -> dropped
)


def _moonshot_fetch():
    return MOONSHOT_CSV


def _chembl_fetch(pages=6):
    return [
        {"molecule_chembl_id": "CHEMBL_A", "canonical_smiles": "c1ccccc1",
         "pchembl_value": "7.2", "standard_value": "60", "activity_comment": None},
        {"molecule_chembl_id": "CHEMBL_I", "canonical_smiles": "c1ccncc1",
         "pchembl_value": "4.0", "standard_value": "30000", "activity_comment": None},
        {"molecule_chembl_id": "CHEMBL_CMT", "canonical_smiles": "c1ccccc1O",
         "pchembl_value": None, "standard_value": None, "activity_comment": "Not Active"},
        # Duplicate of the Moonshot active CCO -> must be deduped (Moonshot kept).
        {"molecule_chembl_id": "CHEMBL_DUP", "canonical_smiles": "OCC",
         "pchembl_value": "6.5", "standard_value": "200", "activity_comment": None},
    ]


def test_parse_moonshot_labels():
    recs = parse_moonshot(MOONSHOT_CSV)
    by_id = {r["id"]: r for r in recs}
    assert by_id["ACT-1"]["label"] == "active"
    assert by_id["INA-1"]["label"] == "inactive"
    assert by_id["INA-2"]["label"] == "inactive"
    assert by_id["INT-1"]["label"] == "intermediate"
    assert "DROP-1" not in by_id  # no measurement -> dropped, not fabricated
    assert all(r["source"] == "moonshot" for r in recs)


def test_parse_chembl_labels():
    recs = parse_chembl(_chembl_fetch())
    by_id = {r["id"]: r for r in recs}
    assert by_id["CHEMBL_A"]["label"] == "active"
    assert by_id["CHEMBL_I"]["label"] == "inactive"
    assert by_id["CHEMBL_CMT"]["label"] == "inactive"  # explicit "Not Active" comment
    assert all(r["source"] == "chembl" for r in recs)


def test_canonical_smiles_dedup_key():
    assert canonical_smiles("OCC") == canonical_smiles("CCO")


def test_assemble_merges_dedups_and_counts():
    ds = assemble(moonshot_fetch=_moonshot_fetch, chembl_fetch=_chembl_fetch)
    # CCO appears in both Moonshot (ACT-1) and ChEMBL (CHEMBL_DUP) -> deduped, Moonshot kept.
    active_ids = {r["id"] for r in ds["actives"]}
    assert "ACT-1" in active_ids
    assert "CHEMBL_DUP" not in active_ids  # Moonshot preferred over the ChEMBL duplicate
    assert ds["n_actives"] >= 2 and ds["n_inactives"] >= 2
    assert ds["actives_by_source"].get("moonshot", 0) >= 1
    assert ds["actives_by_source"].get("chembl", 0) >= 1  # CHEMBL_A is a unique active
    # Provenance records both sources.
    assert {s["source"] for s in ds["sources"]} == {"moonshot", "chembl"}


def test_assemble_degrades_when_moonshot_fails():
    def boom():
        raise RuntimeError("moonshot down")
    ds = assemble(moonshot_fetch=boom, chembl_fetch=_chembl_fetch)
    # No crash; ChEMBL-primary, and the failure is recorded in provenance.
    moon = next(s for s in ds["sources"] if s["source"] == "moonshot")
    assert "error" in moon
    assert ds["n_actives"] >= 1  # ChEMBL still contributes
