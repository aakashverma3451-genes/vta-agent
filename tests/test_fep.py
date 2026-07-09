"""FEP / ABFE seam tests — hermetic, no alchemical MD runner."""
from __future__ import annotations

import pytest

import vta.nodes.fep as fep
from vta.state import new_state


def _state_with_candidates(tmp_path):
    pdb = tmp_path / "pb1.pdb"
    pdb.write_text("END\n")
    st = new_state("x", "run1")
    st["structures"] = {"PB1": {"pdb_path": str(pdb)}}
    st["md_validated_leads"] = [
        {"ligand": "A", "protein": "PB1", "smiles": "CCO", "md_score": 0.8},
        {"ligand": "B", "protein": "PB1", "smiles": "CCN", "md_score": 0.7},
        {"ligand": "C", "protein": "PB1", "smiles": "CCC", "md_score": 0.6},
        {"ligand": "D", "protein": "PB1", "smiles": "CCCl", "md_score": 0.5},
    ]
    return st


def test_fep_skips_when_runner_absent(monkeypatch, tmp_path):
    monkeypatch.setattr(fep, "_fep_runner", lambda: None)
    out = fep.fep_node(_state_with_candidates(tmp_path))
    assert len(out["fep_results"]) == 3
    assert all(r["status"] == "skipped" for r in out["fep_results"].values())
    assert all(c["fep_badge"] == "FEP not run" for c in out["fep_validated_leads"])
    assert any("[skip] FEP" in line for line in out["audit_trail"])


def test_fep_runs_and_sorts_by_delta_g(monkeypatch, tmp_path):
    monkeypatch.setattr(fep, "_fep_runner", lambda: "/usr/local/bin/openfe-abfe")

    values = {"A": -7.0, "B": -9.5, "C": -6.2}

    def fake_run(_cmd, _pdb, smiles, _out_dir):
        ligand = {"CCO": "A", "CCN": "B", "CCC": "C"}[smiles]
        return {"status": "completed", "delta_g": values[ligand], "error": 0.5}

    monkeypatch.setattr(fep, "_run_abfe", fake_run)
    out = fep.fep_node(_state_with_candidates(tmp_path))
    assert [c["ligand"] for c in out["fep_validated_leads"]] == ["B", "A", "C"]
    assert out["fep_validated_leads"][0]["fep_delta_g"] == -9.5
    assert any("FEP: 3/3 candidates completed" in line for line in out["audit_trail"])


def test_fep_no_candidates():
    out = fep.fep_node(new_state("x", "x"))
    assert out["fep_results"] == {}
    assert out["fep_validated_leads"] == []


def test_graph_compiles_with_fep():
    pytest.importorskip("langgraph", reason="requires langgraph")
    from vta.graph import build_app

    app = build_app(include_fep=True)
    assert app is not None
