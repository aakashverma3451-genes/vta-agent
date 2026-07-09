"""Phase 9 actives and bootstrap contracts."""
from __future__ import annotations

from scripts.fetch_phase9_actives import collapse_actives
from scripts.build_hcv_ns5b_ni_scope import build_scope
from scripts.stratify_hcv_ns5b_actives import classify


def test_collapse_actives_keeps_best_per_molecule():
    rows = [
        {
            "molecule_chembl_id": "CHEMBL1",
            "canonical_smiles": "CCO",
            "molecule_pref_name": "A",
            "standard_value": "100",
            "standard_type": "IC50",
            "standard_units": "nM",
            "standard_relation": "=",
            "pchembl_value": "7.0",
            "assay_chembl_id": "ASSAY1",
            "document_chembl_id": "DOC1",
            "target_chembl_id": "TARGET1",
        },
        {
            "molecule_chembl_id": "CHEMBL1",
            "canonical_smiles": "CCO",
            "molecule_pref_name": "A",
            "standard_value": "10",
            "standard_type": "IC50",
            "standard_units": "nM",
            "standard_relation": "=",
            "pchembl_value": "8.0",
            "assay_chembl_id": "ASSAY2",
            "document_chembl_id": "DOC2",
            "target_chembl_id": "TARGET1",
        },
    ]
    out = collapse_actives(rows)
    assert len(out) == 1
    assert out[0]["best_pchembl"] == 8.0
    assert out[0]["assay_chembl_id"] == "ASSAY2"


def test_collapse_actives_filters_weak_or_unqualified_rows():
    rows = [
        {"molecule_chembl_id": "A", "standard_value": "20000", "pchembl_value": "4", "standard_relation": "="},
        {"molecule_chembl_id": "B", "standard_value": "100", "pchembl_value": "7", "standard_relation": ">"},
    ]
    assert collapse_actives(rows) == []


def test_hcv_stratifier_excludes_known_off_target():
    out = classify({"pref_name": "Elbasvir", "canonical_smiles": "CC", "target_chembl_id": "CHEMBL5375"}, [])
    assert out["mechanism_class"] == "off_target_hcv_drug"
    assert out["benchmark_scope"] == "exclude"


def test_hcv_stratifier_labels_filibuvir_assay_as_nni():
    out = classify({
        "pref_name": "CHEMBL230744",
        "canonical_smiles": "CC",
        "target_chembl_id": "CHEMBL5375",
        "assay_chembl_id": "CHEMBL5348065",
    }, [])
    assert out["mechanism_class"] == "NNI"
    assert out["benchmark_scope"] == "allosteric"


def test_hcv_ni_scope_artifact_has_powered_active_count():
    out = build_scope()
    assert out["mechanism_class"] == "NI"
    assert out["n_actives"] >= 20
