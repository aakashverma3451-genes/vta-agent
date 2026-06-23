"""MD pipeline tests — fully hermetic (no OpenMM, no MDAnalysis, no GPU)."""
from __future__ import annotations

import vta.nodes.md_analyze as md_analyze
import vta.nodes.md_simulate as md_simulate
from vta.nodes.md_rerank import md_rerank_node
from vta.nodes.md_select import md_select_node
from vta.state import new_state


# ── helpers ────────────────────────────────────────────────────────────────────
def _leads(n=5):
    return [
        {"ligand": f"Drug{i}", "protein": "PB1", "score": round(0.9 - i * 0.1, 2),
         "smiles": "C1=NC2=NC=NC2=N1", "positive_control": False}
        for i in range(n)
    ]


def _st_with_leads(n=5):
    st = new_state("x", "x")
    st["lead_candidates"] = _leads(n)
    return st


# ── md_select ──────────────────────────────────────────────────────────────────
def test_select_picks_top5():
    st = md_select_node(_st_with_leads(8))
    assert len(st["md_candidates"]) == 5
    assert st["md_candidates"][0]["ligand"] == "Drug0"
    assert any("MD-select" in l for l in st["audit_trail"])


def test_select_noop_when_no_leads():
    st = md_select_node(new_state("x", "x"))
    assert st["md_candidates"] == []
    assert any("MD-select: no leads" in l for l in st["audit_trail"])


# ── md_simulate ────────────────────────────────────────────────────────────────
def test_simulate_skips_gracefully_no_openmm(monkeypatch):
    monkeypatch.setattr(md_simulate, "_check_openmm", lambda: False)
    st = new_state("x", "x")
    st["md_candidates"] = _leads(2)
    out = md_simulate.md_simulate_node(st)
    assert all(v["status"] == "skipped" for v in out["md_results"].values())
    assert any("[skip] MD-simulate" in l for l in out["audit_trail"])
    assert "skipped" in out["versions"]["md"]


def test_simulate_noop_when_no_candidates(monkeypatch):
    monkeypatch.setattr(md_simulate, "_check_openmm", lambda: True)
    st = new_state("x", "x")
    st["md_candidates"] = []
    out = md_simulate.md_simulate_node(st)
    assert out["md_results"] == {}


# ── md_analyze ─────────────────────────────────────────────────────────────────
def test_analyze_skips_gracefully_no_mdanalysis(monkeypatch):
    monkeypatch.setattr(md_analyze, "_check_mdanalysis", lambda: False)
    st = new_state("x", "x")
    st["md_results"] = {"Drug0": {"status": "completed",
                                   "trajectory": "/t.dcd", "final_pdb": "/f.pdb"}}
    out = md_analyze.md_analyze_node(st)
    assert all(v["status"] == "skipped" for v in out["md_analysis"].values())
    assert any("[skip] MD-analyze" in l for l in out["audit_trail"])


def test_analyze_skips_failed_simulations(monkeypatch):
    monkeypatch.setattr(md_analyze, "_check_mdanalysis", lambda: True)
    st = new_state("x", "x")
    st["md_results"] = {"Drug0": {"status": "failed", "reason": "OOM"}}
    out = md_analyze.md_analyze_node(st)
    assert out["md_analysis"]["Drug0"]["status"] == "skipped"


def test_analyze_noop_empty_results():
    st = new_state("x", "x")
    st["md_results"] = {}
    out = md_analyze.md_analyze_node(st)
    assert out["md_analysis"] == {}


# ── md_rerank ──────────────────────────────────────────────────────────────────
def _rerank_state(verdicts: dict):
    st = new_state("x", "x")
    st["md_candidates"] = [
        {"ligand": lid, "score": 0.6, "protein": "PB1"}
        for lid in verdicts
    ]
    st["md_analysis"] = {
        lid: {"md_verdict": v} for lid, v in verdicts.items()
    }
    return st


def test_rerank_stable_gets_boost():
    st = md_rerank_node(_rerank_state({"Riba": "STABLE"}))
    r = st["md_validated_leads"][0]
    assert r["md_score"] > r["score"]
    assert r["md_badge"] == "MD-STABLE"


def test_rerank_unstable_gets_penalty():
    st = md_rerank_node(_rerank_state({"Lopi": "UNSTABLE"}))
    r = st["md_validated_leads"][0]
    assert r["md_score"] < r["score"]
    assert r["md_badge"] == "MD-UNSTABLE"


def test_rerank_sorts_by_md_score():
    st = md_rerank_node(_rerank_state(
        {"A": "UNSTABLE", "B": "STABLE", "C": "MODERATE"}))
    badges = [r["md_badge"] for r in st["md_validated_leads"]]
    assert badges[0] == "MD-STABLE"
    assert badges[-1] == "MD-UNSTABLE"


def test_rerank_not_run_keeps_original_score():
    st = md_rerank_node(_rerank_state({"X": "not_run"}))
    r = st["md_validated_leads"][0]
    assert r["md_score"] == r["score"]
    assert r["md_badge"] == "MD not run"


# ── MM-GBSA binding energy ───────────────────────────────────────────────────
# A trimmed FINAL_RESULTS_MMPBSA.dat as written by MMPBSA.py / gmx_MMPBSA.
_SAMPLE_DAT = """\
Differences (Complex - Receptor - Ligand):
Energy Component        Average      Std. Dev.   Std. Err. of Mean
-------------------------------------------------------------------------------
VDWAALS                -41.2345        2.1100        0.0700
EEL                    -15.6000        3.0000        0.0900
EGB                     30.1000        1.5000        0.0500
DELTA TOTAL            -28.4500        3.1200        0.1000
"""


def test_parse_mmgbsa_dat_extracts_delta_total():
    parsed = md_analyze._parse_mmgbsa_dat(_SAMPLE_DAT)
    assert parsed["mean_binding_energy"] == -28.45
    assert parsed["std_binding_energy"] == 3.12


def test_parse_mmgbsa_dat_handles_unicode_delta():
    # gmx_MMPBSA renders the difference row as 'ΔTOTAL'.
    txt = "ΔTOTAL            -30.10        2.50        0.08\n"
    parsed = md_analyze._parse_mmgbsa_dat(txt)
    assert parsed["mean_binding_energy"] == -30.10
    assert parsed["std_binding_energy"] == 2.50


def test_parse_mmgbsa_dat_returns_none_when_absent():
    assert md_analyze._parse_mmgbsa_dat("no binding line here\n") is None


def test_mmgbsa_skips_when_engine_absent(monkeypatch):
    monkeypatch.setattr(md_analyze, "_mmgbsa_bin", lambda: None)
    out = md_analyze._compute_mmgbsa(
        {"parm": "/x/c.prmtop", "trajectory": "/x/t.dcd"}, "/x")
    assert out["mean_binding_energy"] is None
    assert "[skip] MM-GBSA" in out["note"]


def test_mmgbsa_notes_missing_topology(monkeypatch):
    # Engine present but md_simulate hasn't emitted an Amber prmtop yet.
    monkeypatch.setattr(md_analyze, "_mmgbsa_bin", lambda: "gmx_MMPBSA")
    out = md_analyze._compute_mmgbsa({"trajectory": "/x/t.dcd"}, "/x")
    assert out["mean_binding_energy"] is None
    assert "prmtop" in out["note"]


def test_mmgbsa_computes_when_engine_available(monkeypatch):
    monkeypatch.setattr(md_analyze, "_mmgbsa_bin", lambda: "gmx_MMPBSA")
    monkeypatch.setattr(md_analyze, "_run_mmgbsa_tool",
                        lambda b, p, t, d: _SAMPLE_DAT)
    out = md_analyze._compute_mmgbsa(
        {"parm": "/x/c.prmtop", "trajectory": "/x/t.dcd"}, "/x")
    assert out["mean_binding_energy"] == -28.45
    assert out["std_binding_energy"] == 3.12
    assert out["method"] == "gmx_MMPBSA"


def test_mmgbsa_handles_tool_failure(monkeypatch):
    monkeypatch.setattr(md_analyze, "_mmgbsa_bin", lambda: "MMPBSA.py")
    monkeypatch.setattr(md_analyze, "_run_mmgbsa_tool", lambda b, p, t, d: None)
    out = md_analyze._compute_mmgbsa(
        {"parm": "/x/c.prmtop", "trajectory": "/x/t.dcd"}, "/x")
    assert out["mean_binding_energy"] is None
    assert "no parseable results" in out["note"]


def test_analyze_node_populates_mmgbsa_energy(monkeypatch):
    # Node-level: stub the MDAnalysis-dependent steps + the MM-GBSA engine so the
    # whole md_analyze_node runs hermetically and surfaces a real ΔG.
    monkeypatch.setattr(md_analyze, "_check_mdanalysis", lambda: True)
    monkeypatch.setattr(md_analyze, "_rmsd_analysis",
                        lambda traj, top: {"mean_rmsd": 1.2, "stability": "STABLE"})
    monkeypatch.setattr(md_analyze, "_contacts_analysis",
                        lambda traj, top: {"n_persistent": 3, "top_contacts": []})
    monkeypatch.setattr(md_analyze, "_mmgbsa_bin", lambda: "gmx_MMPBSA")
    monkeypatch.setattr(md_analyze, "_run_mmgbsa_tool",
                        lambda b, p, t, d: _SAMPLE_DAT)

    st = new_state("x", "x")
    st["md_results"] = {"Riba": {"status": "completed", "trajectory": "/x/t.dcd",
                                  "final_pdb": "/x/f.pdb", "parm": "/x/c.prmtop"}}
    out = md_analyze.md_analyze_node(st)
    rec = out["md_analysis"]["Riba"]
    assert rec["md_verdict"] == "STABLE"
    assert rec["mmgbsa"]["mean_binding_energy"] == -28.45
    assert any("MM-GBSA -28.45 kcal/mol" in l for l in out["audit_trail"])


# ── full mock chain ────────────────────────────────────────────────────────────
def test_mock_md_chain_end_to_end(monkeypatch):
    """All 4 MD nodes chain together without crashing when tools are absent."""
    monkeypatch.setattr(md_simulate, "_check_openmm", lambda: False)
    monkeypatch.setattr(md_analyze, "_check_mdanalysis", lambda: False)

    st = _st_with_leads(5)
    st = md_select_node(st)
    st = md_simulate.md_simulate_node(st)
    st = md_analyze.md_analyze_node(st)
    st = md_rerank_node(st)

    assert len(st["md_validated_leads"]) == 5
    assert all(r["md_badge"] == "MD not run" for r in st["md_validated_leads"])
    markers = ["MD-select", "MD-simulate", "MD-analyze", "MD-rerank"]
    trail = "\n".join(st["audit_trail"])
    for m in markers:
        assert m in trail, f"missing audit marker: {m}"
