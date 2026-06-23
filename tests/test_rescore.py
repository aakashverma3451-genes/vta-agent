"""DL-rescore node tests — hermetic (no gnina binary, no GPU).

The real GNINA call is isolated behind `_gnina_bin` / `_run_gnina`, so these tests
monkeypatch those seams and never touch a binary — the same discipline that keeps
the ADMET / Vina tests offline.
"""
from __future__ import annotations

import vta.nodes.rescore as rescore
from vta.state import new_state


def _state(rows):
    st = new_state("x", "x")
    st["docking_results"] = rows
    return st


def test_rescore_annotates_poses_when_gnina_present(monkeypatch, tmp_path):
    pose = tmp_path / "p.pdbqt"; pose.write_text("x")
    rec = tmp_path / "r.pdbqt"; rec.write_text("x")
    rows = [
        {"ligand": "Ribavirin", "dG": -7.1,
         "pose_path": str(pose), "receptor_path": str(rec)},
        {"ligand": "Sofosbuvir", "dG": -8.0,
         "pose_path": str(pose), "receptor_path": str(rec)},
    ]
    monkeypatch.setattr(rescore, "_gnina_bin", lambda: "gnina")
    monkeypatch.setattr(rescore, "_run_gnina",
                        lambda g, r, p: (0.83, 6.2 if "Sofos" not in p else 6.2))

    out = rescore.rescore_node(_state(rows))
    r0 = out["docking_results"][0]
    assert r0["cnn_score"] == 0.83 and r0["cnn_affinity"] == 6.2
    assert out["versions"]["rescore"] == "GNINA (CNN)"
    assert any("DL-rescore: GNINA CNN re-scored 2 poses" in l for l in out["audit_trail"])
    # annotation-only: it must NOT have ranked or dropped anything
    assert len(out["docking_results"]) == 2 and "score" not in r0


def test_rescore_skips_gracefully_when_gnina_absent(monkeypatch):
    rows = [{"ligand": "Ribavirin", "dG": -7.1}]
    monkeypatch.setattr(rescore, "_gnina_bin", lambda: None)

    out = rescore.rescore_node(_state(rows))
    assert "cnn_score" not in out["docking_results"][0]
    assert out["versions"]["rescore"] == "skipped"
    assert any("[skip] DL-rescore: gnina not installed" in l for l in out["audit_trail"])


def test_rescore_skips_records_without_saved_poses(monkeypatch):
    # gnina present but mock-docking records carry no pose_path → nothing to score.
    rows = [{"ligand": "Ribavirin", "dG": -7.1}]
    monkeypatch.setattr(rescore, "_gnina_bin", lambda: "gnina")

    out = rescore.rescore_node(_state(rows))
    assert "cnn_score" not in out["docking_results"][0]
    assert out["versions"]["rescore"] == "skipped (no poses)"
    assert any("no docked poses to score" in l for l in out["audit_trail"])


def test_rescore_noop_when_no_docking_results():
    out = rescore.rescore_node(new_state("x", "x"))
    assert out.get("versions", {}).get("rescore") is None
