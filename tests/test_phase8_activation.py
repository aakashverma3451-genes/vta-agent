"""Phase 8 activation artifact tests."""
from __future__ import annotations

import json

from scripts.phase8_activate import _entries, _verdict
from scripts.fetch_phase8_decoys import _locally_matched, _props


def test_phase8_entries_and_verdict():
    rows = [{"ligand": "A", "score": 1.0, "positive_control": True, "dG": -7, "smiles": "CCO"}]
    assert _entries(rows)[0]["name"] == "A"
    report = {"n_decoys": 9}
    assert "property-matched" in _verdict(report, {"property_matched": False})
    assert "underpowered" in _verdict(report, {"property_matched": True})


def test_local_property_match():
    active = _props("CCO")
    assert _locally_matched(active, "CCN") is True
    assert _locally_matched(active, "CCCCCCCCCCCCCCCCCCCC") is False


def test_phase8_artifacts_exist_after_activation():
    bench = json.load(open("outputs/phase8/matched_benchmark_tilv_pb1.json"))
    assert bench["benchmark"]["n_actives"] == 4
    assert bench["benchmark"]["n_decoys"] >= 1
    locked = json.load(open("outputs/phase8/locked_predictions_tilv_pb1.json"))
    assert locked["status"] == "locked"
