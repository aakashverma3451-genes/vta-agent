"""Phase 2-7 scientific architecture contract tests."""
from __future__ import annotations

import pytest

from vta.data.licenses import validate_registry
from vta.eval.leakage import gate_recalibration_criteria, leakage_audit
from vta.eval.metrics import enrichment_report, log_auc
from vta.eval.splits import scaffold_split, time_split
from vta.experiment.registry import ingest_result, lock_predictions
from vta.kg.schema import validate_edge
from vta.lit.extraction import extract_relation
from vta.nodes.consensus import consensus_node
from vta.nodes.resistance import resistance_node
from vta.nodes.selectivity import selectivity_node
from vta.nodes.target_prioritization import target_prioritization_node
from vta.state import new_state
from vta.structure.qc import qc_structure


def test_phase2_metrics_splits_and_leakage():
    labels = [1, 0, 1, 0, 0]
    assert log_auc(labels) > 0
    rep = enrichment_report([
        {"name": "a", "score": 1.0, "positive_control": True},
        {"name": "d", "score": 0.1, "positive_control": False},
    ])
    assert "log_auc" in rep
    entries = [{"smiles": "CCO", "year": 2020}, {"smiles": "c1ccccc1", "year": 2025}]
    assert set(scaffold_split(entries)) == {"train", "test"}
    assert len(time_split(entries, 2021)["test"]) == 1
    assert leakage_audit(["P12345"], ["p12345"])["promotable"] is False
    assert "logAUC" in gate_recalibration_criteria()["primary_metrics"]


def test_consensus_annotation_only():
    st = new_state("x", "r")
    st["docking_results"] = [
        {"ligand": "A", "dG": -8.0, "cnn_affinity": 6.0, "boltzina_score": 0.8},
        {"ligand": "B", "dG": -6.0, "cnn_affinity": 4.0, "boltzina_score": 0.2},
    ]
    out = consensus_node(st)
    assert out["docking_results"][0]["consensus"]["ranking_active"] is False


def test_selectivity_resistance_and_registry_contracts():
    st = new_state("x", "r")
    st["lead_candidates"] = [{"ligand": "A", "dG": -8.0, "score": 1.0}]
    out = resistance_node(selectivity_node(st))
    assert out["lead_candidates"][0]["selectivity"]["status"] == "skipped_no_counter_target_panel"
    reg = lock_predictions(out["lead_candidates"], "r")
    reg2 = ingest_result(reg, "A", {"EC50": 2, "CC50": 20})
    assert reg2["results"]["A"]["selectivity_index"] == 10
    with pytest.raises(ValueError):
        ingest_result(reg2, "A", {"EC50": 1})


def test_target_prioritization_and_structure_qc():
    st = new_state("x", "r")
    st["candidate_targets"] = [
        {"target": "A", "druggability": 1, "resistance_liability": 0, "essentiality": 1},
        {"target": "B", "druggability": 0, "resistance_liability": 1, "essentiality": 0},
    ]
    out = target_prioritization_node(st)
    assert out["target_prioritization"][0]["target"] == "A"
    assert qc_structure({"pdb_path": None})["status"] == "flagged"


def test_kg_lit_and_license_guards():
    ok, missing = validate_edge({
        "type": "INHIBITS",
        "source": "PMID:1",
        "confidence": 0.9,
        "evidence_span": "Drug A inhibits target B.",
    })
    assert ok, missing
    assert extract_relation({
        "subject": "Drug A",
        "subject_id": "CHEMBL:1",
        "type": "INHIBITS",
        "object": "Target B",
        "object_id": "UniProt:P1",
        "source": "PMID:1",
        "evidence_span": "Drug A did not inhibit target B.",
    }) is None
    assert validate_registry()[0] is True
