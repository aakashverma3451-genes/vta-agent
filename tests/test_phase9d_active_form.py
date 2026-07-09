"""Hermetic tests for Phase 9D/9E pure logic (no docking, no network)."""
from __future__ import annotations

import json

from scripts.run_phase9d_active_form import _largest_fragment, gdd_center, is_triphosphate
import scripts.freeze_benchmark as freeze


def test_is_triphosphate_detects_motif():
    tri = "C[C@@]1(O)[C@H](O)[C@@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O[C@H]1n1ccc(=O)[nH]c1=O"
    parent = "CC(C)OC(=O)[C@H](C)N[P@](=O)(OC[C@H]1O[C@@H](n2ccc(=O)[nH]c2=O)[C@](C)(F)[C@@H]1O)Oc1ccccc1"
    assert is_triphosphate(tri) is True
    assert is_triphosphate(parent) is False


def test_largest_fragment_strips_counterion():
    salted = ("C#C[C@@]1(O)[C@H](O)[C@@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O[C@H]1"
              "n1ccc(=O)[nH]c1=O.CCN(CC)CC.CCN(CC)CC")
    out = _largest_fragment(salted)
    assert "CCN(CC)CC" not in out
    assert is_triphosphate(out)


def test_gdd_center_is_aspartate_centroid():
    # Two of the three GDD aspartate CA atoms; centroid is their average.
    pdb = (
        "ATOM      1  CA  ASP A 220      0.000   0.000   0.000  1.00 0.00           C\n"
        "ATOM      2  CA  ASP A 318      6.000   0.000   0.000  1.00 0.00           C\n"
        "ATOM      3  CA  ASP A 319      0.000   6.000   0.000  1.00 0.00           C\n"
        "ATOM      4  CA  ALA A 100     99.000  99.000  99.000  1.00 0.00           C\n"  # ignored
        "END\n"
    )
    c = gdd_center(pdb, chain="A")
    assert c == [2.0, 2.0, 0.0]


def test_freeze_bundles_targets_and_hashes(tmp_path, monkeypatch):
    # Minimal fake artifacts for the freeze builder.
    mpro = {
        "structure": "7L11 chain A", "control_set": "measured inactives",
        "binding_mode": "non_covalent", "real_vina": True, "receptor_prep": "openbabel_fallback",
        "benchmark": {"n": 4}, "verdict": "powered benchmark",
        "selection": {"seed": 17},
        "bootstrap": {"bootstrap": {"bedroc": {"median": 0.68, "ci95": [0.36, 0.89]},
                                     "roc_auc": {"median": 0.58, "ci95": [0.47, 0.69]}}},
        "rows": [
            {"ligand_id": "A1", "positive_control": True},
            {"ligand_id": "A2", "positive_control": True},
            {"ligand_id": "D1", "positive_control": False},
        ],
    }
    (tmp_path / "mpro.json").write_text(json.dumps(mpro))
    monkeypatch.setattr(freeze, "MPRO", tmp_path / "mpro.json")
    monkeypatch.setattr(freeze, "TILV", tmp_path / "absent_tilv.json")
    monkeypatch.setattr(freeze, "HCV_QUALITY", tmp_path / "absent_q.json")
    monkeypatch.setattr(freeze, "HCV_PROV", tmp_path / "absent_p.json")
    monkeypatch.setattr(freeze, "P9D", tmp_path / "absent_9d.json")

    payload = freeze.build()
    mpro_t = next(t for t in payload["targets"] if "Mpro" in t["target"])
    assert mpro_t["splits"]["n_actives"] == 2
    assert mpro_t["splits"]["n_inactives"] == 1
    assert mpro_t["splits"]["active_ids"] == ["A1", "A2"]
    assert mpro_t["metrics_ci"]["ROC_AUC"]["ci95"] == [0.47, 0.69]
    assert len(payload["content_hash"]) == 16
    assert payload["gate_decision"].startswith("DO NOT PROMOTE")


def test_freeze_hash_is_deterministic(tmp_path, monkeypatch):
    for attr in ("MPRO", "TILV", "HCV_QUALITY", "HCV_PROV", "P9D"):
        monkeypatch.setattr(freeze, attr, tmp_path / f"absent_{attr}.json")
    h1 = freeze.build()["content_hash"]
    h2 = freeze.build()["content_hash"]
    assert h1 == h2
