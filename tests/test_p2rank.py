"""P2Rank consensus tests — hermetic (no Java / P2Rank binary)."""
from __future__ import annotations

import vta.nodes.pockets as pockets
from vta.state import new_state


def _base_state():
    st = new_state("x", "x")
    st["structures"] = {"PA": {"pdb_path": "/fake/PA.pdb"}}
    return st


def _fake_fpocket(name, pdb, fpocket):
    return [
        {"id": 1, "center": [10.0, 20.0, 30.0], "volume": 500.0,
         "druggability": 0.75, "score": 0.8, "conservation": 0.5},
        {"id": 2, "center": [50.0, 50.0, 50.0], "volume": 300.0,
         "druggability": 0.50, "score": 0.5, "conservation": 0.5},
    ]


def test_p2rank_adds_consensus_when_nearby(monkeypatch, tmp_path):
    monkeypatch.setattr(pockets, "_fpocket_bin", lambda: "fpocket")
    monkeypatch.setattr(pockets, "_p2rank_bin", lambda: "prank")
    monkeypatch.setattr(pockets, "_pockets_fpocket", _fake_fpocket)
    # P2Rank finds a pocket right next to pocket 1 (dist ≈ 0 Å) — consensus.
    monkeypatch.setattr(pockets, "_pockets_p2rank",
                        lambda pdb, p2rank, wd: [
                            {"rank": 1, "score": 0.9, "probability": 0.85,
                             "center": [10.1, 20.1, 30.1]},
                        ])
    st = pockets.pockets_node(_base_state())
    pa = st["pockets"]["PA"]
    assert pa[0]["consensus"] is True      # pocket 1 agrees
    assert pa[1]["consensus"] is False     # pocket 2 too far away
    assert "p2rank" in pa[0]["detectors"]
    assert any("P2Rank[PA]" in l for l in st["audit_trail"])
    assert "p2rank" in st["versions"]["pockets"]


def test_p2rank_absent_pockets_get_consensus_false(monkeypatch):
    monkeypatch.setattr(pockets, "_fpocket_bin", lambda: "fpocket")
    monkeypatch.setattr(pockets, "_p2rank_bin", lambda: None)
    monkeypatch.setattr(pockets, "_pockets_fpocket", _fake_fpocket)
    st = pockets.pockets_node(_base_state())
    pa = st["pockets"]["PA"]
    for p in pa:
        assert p["consensus"] is False
    assert "p2rank" not in st["versions"].get("pockets", "")


def test_p2rank_failure_degrades_gracefully(monkeypatch):
    monkeypatch.setattr(pockets, "_fpocket_bin", lambda: "fpocket")
    monkeypatch.setattr(pockets, "_p2rank_bin", lambda: "prank")
    monkeypatch.setattr(pockets, "_pockets_fpocket", _fake_fpocket)
    monkeypatch.setattr(pockets, "_pockets_p2rank",
                        lambda *a, **kw: (_ for _ in ()).throw(RuntimeError("P2Rank crash")))
    st = pockets.pockets_node(_base_state())
    pa = st["pockets"]["PA"]
    for p in pa:
        assert p["consensus"] is False
    assert any("WARNING failed" in l for l in st["audit_trail"])


def test_experimental_site_has_consensus_true(monkeypatch):
    monkeypatch.setattr(pockets, "_fpocket_bin", lambda: None)
    monkeypatch.setattr(pockets, "_p2rank_bin", lambda: None)
    st = new_state("x", "x")
    st["structures"] = {"PB1": {"pdb_path": "/fake/PB1.pdb"}}
    out = pockets.pockets_node(st)
    assert out["pockets"]["PB1"][0]["consensus"] is True
    assert out["pockets"]["PB1"][0]["detectors"] == ["experimental"]
