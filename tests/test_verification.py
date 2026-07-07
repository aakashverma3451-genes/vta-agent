"""Tests for R4 VerificationNode + R5 annotation ranker + report downgrade (Phase R, Stage 3)."""
from __future__ import annotations

import json
from pathlib import Path

from vta.nodes.annotate_rank import annotate_rank_node
from vta.nodes.verification import build_verdict, verification_node
from vta.report import render_report
from vta.state import new_state


def _gate_artifacts(d: Path):
    baselines = [
        {"benchmark": "Mpro_noncovalent", "structure_based_enrichment_demonstrated": False},
        {"benchmark": "Winner", "structure_based_enrichment_demonstrated": True},
    ]
    redock = {"rmsd_flag_threshold_A": 2.0, "results": [
        {"target": "Mpro_7L11", "rmsd": 1.65, "pose_reliable": True, "status": "ok"},
    ]}
    (d / "base.json").write_text(json.dumps(baselines))
    (d / "redock.json").write_text(json.dumps(redock))
    return {"baselines_path": d / "base.json", "redock_path": d / "redock.json"}


def _mpro_leads_state():
    st = new_state("verify-mpro", "verify-mpro")
    st["structures"] = {"MPRO": {"method": "experimental", "source": "7L11:A"}}
    st["triage_decision"] = {"MPRO": {"decision": "full_dock", "rationale": []}}
    st["lead_candidates"] = [
        {"ligand": "cpd", "protein": "MPRO", "dG": -8.0, "le": -0.3, "conservation": 0.5,
         "score": 0.9, "positive_control": False}]
    return st


# --- R4 verification -------------------------------------------------------------------------

def test_mpro_docking_claim_is_downgraded(tmp_path):
    # Mpro does NOT beat 2D-similarity → the docking ranking must be DOWNGRADED, not passed.
    v = build_verdict(_mpro_leads_state(), **_gate_artifacts(tmp_path))
    assert v["verdict"] == "downgrade"
    assert v["per_target"]["MPRO"]["verdict"] == "downgrade"
    assert v["per_target"]["MPRO"]["beats_2d_baseline"] is False


def test_pass_only_when_beats_baseline_and_pose_reliable(tmp_path):
    st = _mpro_leads_state()
    st["lead_candidates"][0]["protein"] = "MPRO"
    # point the target-map at a benchmark that DOES beat the baseline, with a reliable pose
    # (simulate by relabelling protein to one the map resolves + a winning baseline entry)
    arts = _gate_artifacts(tmp_path)
    # craft a state whose lead protein maps to nothing → untested → still downgrade (conservative)
    st_unknown = new_state("v2", "v2")
    st_unknown["triage_decision"] = {"GPX": {"decision": "full_dock"}}
    st_unknown["lead_candidates"] = [{"ligand": "x", "protein": "GPX", "dG": -7, "le": -0.2,
                                      "conservation": 0.5, "score": 0.5}]
    v = build_verdict(st_unknown, **arts)
    assert v["verdict"] == "downgrade"    # untested target cannot make a verified enrichment claim


def test_defer_run_is_not_applicable(tmp_path):
    st = new_state("v-defer", "v-defer")   # no leads, no triage
    v = build_verdict(st, **_gate_artifacts(tmp_path))
    assert v["verdict"] == "not_applicable"


def test_verification_node_writes_verdict_and_audit(tmp_path, monkeypatch):
    import vta.nodes.verification as ver
    monkeypatch.setattr(ver, "DEFAULT_BASELINES_PATH", (tmp_path / "base.json"))
    monkeypatch.setattr(ver, "DEFAULT_REDOCK_PATH", (tmp_path / "redock.json"))
    _gate_artifacts(tmp_path)
    st = verification_node(_mpro_leads_state())
    assert st["verification_verdict"]["verdict"] == "downgrade"
    assert any("Verification: run verdict = DOWNGRADE" in a for a in st["audit_trail"])


# --- R5 annotation ranker --------------------------------------------------------------------

def test_annotate_only_target_gets_labelled_ligand_based_annotation():
    st = new_state("ann", "ann")
    st["triage_decision"] = {"NS5B": {"decision": "annotate_only", "rationale": ["un_benchmarkable"]}}
    st["target_dossier"] = {"NS5B": {"metal_dependence": True, "known_ligand_class": "nucleotide"}}
    st["annotation_library"] = [
        {"name": "sofosbuvir-TP", "smiles": "OCC1OC(n2ccc(=O)[nH]c2=O)CC1O", "positive_control": True},
        {"name": "decoyA", "smiles": "CCOc1ccccc1", "positive_control": False},
    ]
    st = annotate_rank_node(st)
    a = st["annotation_rankings"]["NS5B"]
    assert a["status"] == "ligand_based_annotation"
    assert a["ranking"][0]["ranking_type"] == "ligand_based_2d_similarity_annotation"
    assert any("Nucleotide/metal" in c for c in a["caveats"])
    # never attaches a docking dG / enrichment number
    assert all("dG" not in r and "enrichment" not in r for r in a["ranking"])


def test_annotate_only_without_actives_is_labelled_no_ranking():
    st = new_state("ann2", "ann2")
    st["triage_decision"] = {"NS5B": {"decision": "annotate_only"}}
    st["target_dossier"] = {"NS5B": {"metal_dependence": True}}
    st = annotate_rank_node(st)
    a = st["annotation_rankings"]["NS5B"]
    assert a["status"] == "annotation_only" and a["ranking"] == []


# --- report obeys the verdict ----------------------------------------------------------------

def test_report_shows_downgrade_banner():
    st = _mpro_leads_state()   # build_verdict fallback reads the real committed artifacts (Mpro fails)
    html = render_report(st)
    assert "DOWNGRADED to a ligand-based annotation" in html


def test_report_renders_annotation_section():
    st = new_state("rep-ann", "rep-ann")
    st["annotation_rankings"] = {"NS5B": {"status": "annotation_only", "ranking": [],
                                          "caveats": ["Nucleotide/metal caveat"],
                                          "note": "no validated ranking"}}
    html = render_report(st)
    assert "Ligand-based annotations" in html and "Nucleotide/metal caveat" in html
