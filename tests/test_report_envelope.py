"""Tests for the honesty envelope (D0.5) — the non-removable output contract.

Hermetic: writes fixture artifacts to a tmp dir and points the builder at them, so no test
depends on the live committed outputs/ (which can change). Also asserts the STRUCTURAL
guarantee — that render_report always emits the disclaimer + pinned benchmark, for a full run,
a deferred run, and an empty run — which is the whole point of D0.5.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from vta.nodes.rank import WEIGHTS as LIVE_RANK_WEIGHTS
from vta.report import render_report
from vta.report_envelope import DISCLAIMER, build_envelope
from vta.state import new_state


# --- fixtures: minimal, well-formed committed-artifact stand-ins ----------------------------

def _write_artifacts(d: Path) -> dict:
    bench = {
        "artifact": "VTA-Agent frozen validation benchmark (TEST)",
        "frozen_at": "2026-07-03", "content_hash": "deadbeefcafe0001",
        "gate_decision": "DO NOT PROMOTE (annotation-only)",
        "caveats": ["Mpro docking signal is modest (ROC-AUC CI crosses 0.5)."],
        "targets": [
            {"target": "SARS-CoV-2 Mpro (non-covalent)", "grade": "powered (modest signal)",
             "metrics_ci": {"BEDROC": {"median": 0.68, "ci95": [0.36, 0.89]},
                            "ROC_AUC": {"median": 0.58, "ci95": [0.47, 0.69]}}},
            {"target": "HCV NS5B NI (triphosphate)", "grade": "not benchmarkable (recipe)",
             "metrics_ci": None},
        ],
    }
    baselines = [
        {"benchmark": "Mpro_noncovalent",
         "structure_based_enrichment_demonstrated": False,
         "wi3_plain_language": "Vina does NOT beat 2D-similarity on BEDROC (paired test)."},
    ]
    redock = {"rmsd_flag_threshold_A": 2.0, "results": [
        {"target": "Mpro_7L11", "rmsd": 1.654, "pose_reliable": True, "status": "ok"},
        {"target": "HCV_NS5B_2XI3", "rmsd": None, "pose_reliable": None,
         "status": "RMSD computation failed"},
    ]}
    (d / "bench.json").write_text(json.dumps(bench))
    (d / "base.json").write_text(json.dumps(baselines))
    (d / "redock.json").write_text(json.dumps(redock))
    return {"benchmark_path": d / "bench.json", "baselines_path": d / "base.json",
            "redock_path": d / "redock.json"}


def _mpro_state() -> dict:
    st = new_state("env-mpro", "env-mpro")
    st["structures"] = {"MPRO": {"method": "experimental", "source": "7L11:A"}}
    st["lead_candidates"] = [
        {"ligand": "MAT-POS-x", "protein": "MPRO", "dG": -8.1, "le": -0.3,
         "conservation": 0.5, "score": 0.9, "positive_control": True},
    ]
    return st


# --- builder behaviour ----------------------------------------------------------------------

def test_envelope_pins_benchmark_and_matches_target(tmp_path):
    env = build_envelope(_mpro_state(), **_write_artifacts(tmp_path))
    assert env["disclaimer"] == DISCLAIMER
    assert env["benchmark_pin"]["content_hash"] == "deadbeefcafe0001"
    basis = env["ranking_basis"]
    assert len(basis) == 1 and basis[0]["benchmark_target"] == "SARS-CoV-2 Mpro (non-covalent)"
    # trivial-baseline verdict is carried through, honestly negative
    assert basis[0]["trivial_baseline"]["beats_trivial_2d_baseline"] is False
    # pose reliability from redock artifact
    assert basis[0]["pose_reliability"]["status"] == "pose-reliable"
    assert basis[0]["pose_reliability"]["rmsd_A"] == 1.654


def test_envelope_weights_track_live_ranker(tmp_path):
    env = build_envelope(_mpro_state(), **_write_artifacts(tmp_path))
    # the envelope must report the ACTUAL in-force weights, never a hard-coded copy
    assert env["ranking_weights"] == dict(LIVE_RANK_WEIGHTS)
    assert env["ranking_weights"]["dG"] == 1.0 and env["ranking_weights"]["le"] == 0.0


def test_envelope_flags_out_of_validated_domain(tmp_path):
    st = new_state("env-ood", "env-ood")
    st["structures"] = {"GP120": {"method": "esmfold"}}
    st["lead_candidates"] = [{"ligand": "cpd", "protein": "GP120", "dG": -7.0, "le": -0.3,
                              "conservation": 0.5, "score": 0.5}]
    env = build_envelope(st, **_write_artifacts(tmp_path))
    assert "GP120" in env["out_of_validated_domain"]
    assert env["ranking_basis"] == []


def test_envelope_carries_nucleotide_caveat(tmp_path):
    st = new_state("env-hcv", "env-hcv")
    st["lead_candidates"] = [{"ligand": "sofosbuvir", "protein": "NS5B", "dG": -7.5,
                              "le": -0.2, "conservation": 0.5, "score": 0.4}]
    env = build_envelope(st, **_write_artifacts(tmp_path))
    assert any("NUCLEOTIDE" in c for c in env["scoring_caveats"])


def test_envelope_degrades_when_artifacts_missing(tmp_path):
    # every artifact absent → labelled "not available", but the disclaimer still stands
    env = build_envelope(_mpro_state(),
                         benchmark_path=tmp_path / "nope1.json",
                         baselines_path=tmp_path / "nope2.json",
                         redock_path=tmp_path / "nope3.json")
    assert env["disclaimer"] == DISCLAIMER
    assert env["benchmark_pin"]["available"] is False
    assert env["provenance"]["frozen_benchmark"] == "MISSING"


# --- structural guarantee: no report without the envelope -----------------------------------

@pytest.mark.parametrize("make_state", [
    _mpro_state,
    lambda: {**new_state("env-defer", "env-defer"), "route": "defer",
             "taxon_result": {"species": "Y", "confidence_basis": ""},
             "classification_confidence": 30.0},
    lambda: new_state("env-empty", "env-empty"),
])
def test_report_always_contains_disclaimer(make_state):
    html = render_report(make_state())
    assert "Honesty envelope" in html
    assert "HYPOTHESES" in html            # the non-removable disclaimer text
    assert "clinical, efficacy, or safety" in html


def test_report_node_stores_envelope_on_state(tmp_path, monkeypatch):
    # chdir so report_node's write_report lands under tmp, not the repo's outputs/
    monkeypatch.chdir(tmp_path)
    from vta.nodes.report import report_node
    st = _mpro_state()
    st = report_node(st)
    assert "honesty_envelope" in st and st["honesty_envelope"]["disclaimer"] == DISCLAIMER
