"""Hermetic tests for the WI-7 GNINA rescorer contract (GNINA absent in this env)."""
from __future__ import annotations

from vta.nodes.rescore import rescore_node
from vta.state import new_state


def _state_with_rows():
    st = new_state("rescore_test", "rescore_test")
    st["docking_results"] = [
        {"ligand": "A", "dG": -8.0, "positive_control": True, "score": 0.6},
        {"ligand": "B", "dG": -6.0, "positive_control": False, "score": 0.3},
    ]
    return st


def test_rescore_no_op_when_gnina_absent():
    # No GNINA installed here → clean no-op, ranking untouched, DO NOT PROMOTE logged.
    st = rescore_node(_state_with_rows())
    assert st["versions"]["rescore"] == "skipped"
    assert "DO NOT PROMOTE" in st["rescore_gate"]
    # No cnn_* fields fabricated; docking rows unchanged.
    assert all("cnn_affinity" not in r for r in st["docking_results"])
    assert any("DO NOT PROMOTE" in a for a in st["audit_trail"])


def test_rescore_empty_rows_is_safe():
    st = new_state("empty", "empty")
    st["docking_results"] = []
    out = rescore_node(st)
    assert out["docking_results"] == []  # no crash, nothing invented
