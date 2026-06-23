"""Boltzina DL affinity node tests — hermetic (no GPU, no boltzina package)."""
from __future__ import annotations

import vta.nodes.boltzina as boltzina
from vta.state import new_state


def _state(rows):
    st = new_state("x", "x")
    st["docking_results"] = rows
    st["structures"] = {"PB1": {"pdb_path": "/fake/PB1.pdb"}}
    return st


def _rows():
    return [
        {"ligand": "Ribavirin",  "protein": "PB1", "smiles": "OCC1OC(n2cnc3c(N)ncnc32)C(O)C1O", "dG": -6.8},
        {"ligand": "Sofosbuvir", "protein": "PB1", "smiles": "CC(C)OC(=O)", "dG": -8.5},
    ]


def test_boltzina_annotates_records(monkeypatch):
    monkeypatch.setattr(boltzina, "_boltzina_available", lambda: True)
    monkeypatch.setattr(boltzina, "_boltzina_predict", lambda pdb, smi: 0.82)

    out = boltzina.boltzina_node(_state(_rows()))
    for r in out["docking_results"]:
        assert r["boltzina_score"] == 0.82
    assert out["versions"]["boltzina"].startswith("Boltzina")
    assert any("Boltzina: scored 2 records" in l for l in out["audit_trail"])


def test_boltzina_skips_gracefully_when_unavailable(monkeypatch):
    monkeypatch.setattr(boltzina, "_boltzina_available", lambda: False)

    out = boltzina.boltzina_node(_state(_rows()))
    for r in out["docking_results"]:
        assert "boltzina_score" not in r
    assert out["versions"]["boltzina"] == "skipped"
    assert any("[skip] Boltzina: package not installed" in l for l in out["audit_trail"])


def test_boltzina_skips_records_missing_smiles(monkeypatch):
    monkeypatch.setattr(boltzina, "_boltzina_available", lambda: True)
    monkeypatch.setattr(boltzina, "_boltzina_predict", lambda pdb, smi: 0.7)

    rows = [{"ligand": "NoSmiles", "protein": "PB1", "dG": -7.0}]  # no smiles key
    out = boltzina.boltzina_node(_state(rows))
    assert "boltzina_score" not in out["docking_results"][0]
    assert "skipped" in out["versions"]["boltzina"]


def test_boltzina_noop_when_no_docking_results():
    out = boltzina.boltzina_node(new_state("x", "x"))
    assert out.get("versions", {}).get("boltzina") is None


def test_boltzina_annotation_only_no_ranking_change(monkeypatch):
    """boltzina_score must not appear in or alter the score field."""
    monkeypatch.setattr(boltzina, "_boltzina_available", lambda: True)
    monkeypatch.setattr(boltzina, "_boltzina_predict", lambda pdb, smi: 0.9)

    rows = _rows()
    original_dGs = [r["dG"] for r in rows]
    out = boltzina.boltzina_node(_state(rows))
    assert [r["dG"] for r in out["docking_results"]] == original_dGs
    assert "score" not in out["docking_results"][0]  # rank_node hasn't run yet
