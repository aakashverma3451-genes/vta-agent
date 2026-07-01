"""Hermetic tests for Phase 10 ensemble docking logic (no network, no real Vina)."""
from __future__ import annotations

import json
import os

import pytest

import scripts.run_phase10_ensemble as p10


# ── fixtures ──────────────────────────────────────────────────────────────────
FAKE_PDB_CHAIN_A = (
    "ATOM      1  CA  ALA A   1      10.000  10.000  10.000  1.00 20.00           C\n"
    "ATOM      2  CA  ALA A   2      11.000  10.000  10.000  1.00 20.00           C\n"
    "TER\n"
    "END\n"
)

FAKE_PDB_WITH_DYAD = (
    "ATOM      1  CA  HIS A  41      4.000   6.000   8.000  1.00 20.00           C\n"
    "ATOM      2  CA  CYS A 145     10.000  12.000  14.000  1.00 20.00           C\n"
    "TER\n"
    "END\n"
)


def _phase9_rows():
    return [
        {"ligand_id": "CPDA", "positive_control": True,  "dG": -9.0,  "smiles": "c1ccccc1", "ligand": "CPDA"},
        {"ligand_id": "CPDB", "positive_control": True,  "dG": -8.0,  "smiles": "c1ccccc1", "ligand": "CPDB"},
        {"ligand_id": "CPDC", "positive_control": False, "dG": -5.0,  "smiles": "c1ccccc1", "ligand": "CPDC"},
        {"ligand_id": "CPDD", "positive_control": False, "dG": -4.0,  "smiles": "c1ccccc1", "ligand": "CPDD"},
    ]


def _phase9_json(rows, bench=None):
    return {
        "rows": rows,
        "benchmark": bench or {
            "bedroc": 0.68, "roc_auc": 0.58, "log_auc": 0.20,
            "ef": {"EF1%": 2.0},
        },
        "bootstrap": {
            "bootstrap": {
                "bedroc": {"median": 0.68, "ci95": [0.36, 0.89]},
                "roc_auc": {"median": 0.58, "ci95": [0.47, 0.69]},
            }
        },
    }


# ── merge_ensemble ────────────────────────────────────────────────────────────
def test_merge_takes_best_dg():
    rows = _phase9_rows()
    new_scores = {
        "6Y2E_apo": {"CPDA": -9.5, "CPDB": -7.0, "CPDC": -5.5, "CPDD": -3.0},
        "7K3T_holo": {"CPDA": -8.0, "CPDB": -8.5, "CPDC": -4.8, "CPDD": -4.2},
    }
    merged = p10.merge_ensemble(rows, new_scores)
    assert len(merged) == 4
    # CPDA: primary -9.0, 6Y2E -9.5 → ensemble should pick 6Y2E (-9.5)
    cpda = next(r for r in merged if r["ligand_id"] == "CPDA")
    assert cpda["ensemble_dG"] == pytest.approx(-9.5)
    assert cpda["ensemble_conformer"] == "6Y2E_apo"
    # CPDB: primary -8.0, 6Y2E -7.0, 7K3T -8.5 → primary was best? No: -8.5 < -8.0 → 7K3T
    cpdb = next(r for r in merged if r["ligand_id"] == "CPDB")
    assert cpdb["ensemble_dG"] == pytest.approx(-8.5)
    assert cpdb["ensemble_conformer"] == "7K3T_holo"


def test_merge_falls_back_to_primary_when_new_fail():
    rows = _phase9_rows()
    new_scores = {
        "6Y2E_apo": {"CPDA": None, "CPDB": -7.0, "CPDC": None, "CPDD": -3.0},
    }
    merged = p10.merge_ensemble(rows, new_scores)
    # CPDA new score is None → falls back to primary -9.0
    cpda = next(r for r in merged if r["ligand_id"] == "CPDA")
    assert cpda["ensemble_dG"] == pytest.approx(-9.0)
    assert cpda["ensemble_conformer"] == "7L11_primary"


def test_merge_excludes_compound_when_all_fail():
    rows = [{"ligand_id": "X", "positive_control": True, "dG": None, "smiles": "c1ccccc1", "ligand": "X"}]
    new_scores = {"6Y2E_apo": {"X": None}}
    merged = p10.merge_ensemble(rows, new_scores)
    assert len(merged) == 0


def test_merge_preserves_positive_control_label():
    rows = _phase9_rows()
    new_scores = {"6Y2E_apo": {"CPDA": -9.5, "CPDB": -7.0, "CPDC": -5.5, "CPDD": -3.0}}
    merged = p10.merge_ensemble(rows, new_scores)
    for r in merged:
        expected = next(row["positive_control"] for row in rows if row["ligand_id"] == r["ligand_id"])
        assert r["positive_control"] == expected


# ── _entries_from_rows ────────────────────────────────────────────────────────
def test_entries_excludes_none_dg():
    rows = [
        {"ligand_id": "A", "ensemble_dG": -9.0, "positive_control": True},
        {"ligand_id": "B", "ensemble_dG": None, "positive_control": False},
    ]
    entries = p10._entries_from_rows(rows, "ensemble_dG")
    assert len(entries) == 1
    assert entries[0]["name"] == "A"
    assert entries[0]["score"] == pytest.approx(9.0)


def test_entries_score_is_neg_dg():
    rows = [{"ligand_id": "X", "ensemble_dG": -7.5, "positive_control": True}]
    entries = p10._entries_from_rows(rows)
    assert entries[0]["score"] == pytest.approx(7.5)


# ── _ranking_delta ────────────────────────────────────────────────────────────
def test_ranking_delta_improvement():
    # Phase 9: A ranks 1 (best = most negative dG), B ranks 2
    phase9_rows = [
        {"ligand_id": "A", "positive_control": True,  "dG": -9.0},
        {"ligand_id": "B", "positive_control": True,  "dG": -8.0},
        {"ligand_id": "C", "positive_control": False, "dG": -7.0},
    ]
    # Ensemble: B now has -9.5 (best), A -8.5, C -6.0
    ensemble_rows = [
        {"ligand_id": "A", "ensemble_dG": -8.5, "positive_control": True,  "dG": -9.0},
        {"ligand_id": "B", "ensemble_dG": -9.5, "positive_control": True,  "dG": -8.0},
        {"ligand_id": "C", "ensemble_dG": -6.0, "positive_control": False, "dG": -7.0},
    ]
    delta = p10._ranking_delta(phase9_rows, ensemble_rows)
    a_delta = next(d for d in delta if d["compound_id"] == "A")
    b_delta = next(d for d in delta if d["compound_id"] == "B")
    # A was rank 1, now rank 2 → delta +1, not improved
    assert a_delta["delta_rank"] == 1
    assert a_delta["improved"] is False
    # B was rank 2, now rank 1 → delta -1, improved
    assert b_delta["delta_rank"] == -1
    assert b_delta["improved"] is True


# ── _dyad_center ──────────────────────────────────────────────────────────────
def test_dyad_center_returns_midpoint(tmp_path):
    pdb = tmp_path / "test.pdb"
    pdb.write_text(FAKE_PDB_WITH_DYAD)
    center = p10._dyad_center(str(pdb), chain="A")
    assert center is not None
    assert center == pytest.approx([7.0, 9.0, 11.0])  # (4+10)/2, (6+12)/2, (8+14)/2


def test_dyad_center_returns_none_when_residues_missing(tmp_path):
    pdb = tmp_path / "test.pdb"
    pdb.write_text(FAKE_PDB_CHAIN_A)  # res 1 and 2, not 41/145
    center = p10._dyad_center(str(pdb), chain="A")
    assert center is None


def test_dyad_center_ignores_wrong_chain(tmp_path):
    pdb = tmp_path / "test.pdb"
    # His41 and Cys145 on chain B, not A
    pdb.write_text(
        "ATOM      1  CA  HIS B  41      4.000   6.000   8.000  1.00 20.00           C\n"
        "ATOM      2  CA  CYS B 145     10.000  12.000  14.000  1.00 20.00           C\n"
        "TER\nEND\n"
    )
    center = p10._dyad_center(str(pdb), chain="A")
    assert center is None


# ── fetch_and_prep_conformer (network + openbabel mocked) ─────────────────────
def test_fetch_and_prep_conformer_full_pipeline(tmp_path, monkeypatch):
    """Full pipeline: fetch → extract chain → write PDB → OpenBabel prep → PDBQT."""
    FAKE_FULL_PDB = (
        "ATOM      1  CA  ALA A   1      10.000  10.000  10.000  1.00 20.00           C\n"
        "ATOM      2  CA  ALA B   1      20.000  20.000  20.000  1.00 20.00           C\n"
        "TER\n"
        "END\n"
    )
    monkeypatch.setattr(p10, "_fetch_rcsb", lambda pdb_id: FAKE_FULL_PDB)
    monkeypatch.setattr(p10, "STRUCTURES_DIR", tmp_path)

    fake_pdbqt_content = "REMARK fake receptor\nATOM      1  C   ALA A   1       0.000   0.000   0.000  0.00  0.00           C\n"

    def fake_prep(pdb_path, pdbqt):
        if os.path.exists(pdbqt):
            return pdbqt
        with open(pdbqt, "w") as fh:
            fh.write(fake_pdbqt_content)
        return pdbqt

    monkeypatch.setattr(p10, "_prep_receptor_obabel", fake_prep)

    conf = {"label": "test_apo", "pdb_id": "6Y2E", "chain": "A", "note": "test"}
    receptor, center = p10.fetch_and_prep_conformer(conf)
    assert receptor is not None
    assert receptor.endswith(".pdbqt")
    # PDB should contain only chain A atoms
    pdb_written = (tmp_path / "phase10_6Y2E_A.pdb").read_text()
    assert "ATOM" in pdb_written
    assert "ALA A" in pdb_written
    # center falls back to MPRO_CENTER (no His41/Cys145 in FAKE_PDB_CHAIN_A)
    assert center is not None


def test_fetch_and_prep_conformer_reuses_cache(tmp_path, monkeypatch):
    monkeypatch.setattr(p10, "STRUCTURES_DIR", tmp_path)
    pdbqt_path = tmp_path / "phase10_6Y2E_A_rec.pdbqt"
    pdbqt_path.write_text("REMARK cached\n")

    fetch_called = []
    monkeypatch.setattr(p10, "_fetch_rcsb", lambda x: fetch_called.append(x) or "")
    monkeypatch.setattr(p10, "_prep_receptor_obabel", lambda p, q: str(pdbqt_path))

    # PDB also needs to exist so we don't hit fetch_rcsb
    pdb_path = tmp_path / "phase10_6Y2E_A.pdb"
    pdb_path.write_text(FAKE_PDB_CHAIN_A)

    conf = {"label": "6Y2E_apo", "pdb_id": "6Y2E", "chain": "A", "note": "test"}
    receptor, center = p10.fetch_and_prep_conformer(conf)
    assert receptor == str(pdbqt_path)
    assert fetch_called == []  # network not called


def test_fetch_network_failure_returns_none(tmp_path, monkeypatch):
    monkeypatch.setattr(p10, "STRUCTURES_DIR", tmp_path)
    monkeypatch.setattr(p10, "_fetch_rcsb", lambda x: (_ for _ in ()).throw(OSError("timeout")))
    conf = {"label": "bad", "pdb_id": "XXXX", "chain": "A", "note": ""}
    receptor, center = p10.fetch_and_prep_conformer(conf)
    assert receptor is None
    assert center is None


# ── dock_compounds_against_receptor ──────────────────────────────────────────
def test_dock_compounds_skips_missing_ligand(tmp_path, monkeypatch):
    monkeypatch.setattr(p10, "_LIG_CACHE", str(tmp_path))  # empty dir → no PDBQTs
    monkeypatch.setattr(p10, "STRUCTURES_DIR", tmp_path)
    monkeypatch.setattr(p10, "_run_vina", lambda *a, **k: -8.0)
    compounds = [{"ligand_id": "MISSING", "positive_control": True}]
    scores = p10.dock_compounds_against_receptor(
        compounds, str(tmp_path / "rec.pdbqt"), "test_conf", "vina"
    )
    assert scores["MISSING"] is None


def test_dock_compounds_uses_cached_pdbqt(tmp_path, monkeypatch):
    lig_pdbqt = tmp_path / "LIG1.pdbqt"
    lig_pdbqt.write_text("REMARK fake lig\n")
    monkeypatch.setattr(p10, "_LIG_CACHE", str(tmp_path))
    monkeypatch.setattr(p10, "STRUCTURES_DIR", tmp_path)

    called_with = []
    def fake_vina(vina, receptor, ligand, center, out):
        called_with.append((receptor, ligand))
        return -7.5

    monkeypatch.setattr(p10, "_run_vina", fake_vina)
    compounds = [{"ligand_id": "LIG1", "positive_control": True}]
    scores = p10.dock_compounds_against_receptor(
        compounds, "rec.pdbqt", "6Y2E_apo", "vina"
    )
    assert scores["LIG1"] == pytest.approx(-7.5)
    assert called_with[0][0] == "rec.pdbqt"


# ── ci_overlap logic ──────────────────────────────────────────────────────────
def test_interpret_ci_overlap_message():
    # Two overlapping CIs should mention "stable"
    msg = p10._interpret(
        {"bedroc": 0.68, "roc_auc": 0.58},
        {"bedroc": 0.70, "roc_auc": 0.60},
        [],
        [0.36, 0.89], [0.40, 0.92],
    )
    assert "stable" in msg.lower()


def test_interpret_ci_non_overlap_message():
    # Non-overlapping: Phase 9 CI [0.50, 0.60], ensemble CI [0.70, 0.85]
    msg = p10._interpret(
        {"bedroc": 0.55, "roc_auc": 0.55},
        {"bedroc": 0.75, "roc_auc": 0.75},
        [],
        [0.50, 0.60], [0.70, 0.85],
    )
    assert "not overlap" in msg.lower()


# ── integration: run() with everything mocked ────────────────────────────────
def test_run_integration(tmp_path, monkeypatch):
    """run() smoke test: mocked RCSB + Vina + OpenBabel; verifies JSON output shape."""
    # Write fake phase9 JSON
    rows = _phase9_rows()
    phase9_data = _phase9_json(rows)
    fake_phase9 = tmp_path / "mpro_noncovalent_benchmark.json"
    fake_phase9.write_text(json.dumps(phase9_data))
    monkeypatch.setattr(p10, "PHASE9_JSON", fake_phase9)

    out_json = tmp_path / "ensemble_benchmark.json"
    out_md = tmp_path / "ensemble_benchmark.md"
    monkeypatch.setattr(p10, "OUT_JSON", out_json)
    monkeypatch.setattr(p10, "OUT_MD", out_md)
    monkeypatch.setattr(p10, "OUT_DIR", tmp_path)
    monkeypatch.setattr(p10, "STRUCTURES_DIR", tmp_path)

    # Ligand cache: write fake PDBQTs so dock_compounds finds them
    for c in rows:
        (tmp_path / f"{c['ligand_id']}.pdbqt").write_text("REMARK fake\n")
    monkeypatch.setattr(p10, "_LIG_CACHE", str(tmp_path))

    # Network: fake RCSB
    monkeypatch.setattr(p10, "_fetch_rcsb", lambda x: FAKE_PDB_CHAIN_A)

    # OpenBabel: write a fake PDBQT
    def fake_prep(pdb_path, pdbqt):
        if os.path.exists(pdbqt):
            return pdbqt
        with open(pdbqt, "w") as fh:
            fh.write("REMARK fake receptor\n")
        return pdbqt
    monkeypatch.setattr(p10, "_prep_receptor_obabel", fake_prep)

    # Vina: returns deterministic fake dG
    vina_call_count = [0]
    def fake_vina(vina_bin, receptor, ligand, center, out):
        vina_call_count[0] += 1
        # Return slightly different scores for each conformer
        if "6Y2E" in receptor or "6Y2E" in out:
            return -9.5 if "CPDA" in out else -7.0
        return -8.0
    monkeypatch.setattr(p10, "_run_vina", fake_vina)
    monkeypatch.setattr(p10, "find_tool", lambda *a: "vina")

    result = p10.run(n_resamples=100)

    assert out_json.exists()
    payload = json.loads(out_json.read_text())
    assert "ensemble_benchmark" in payload
    assert "phase9_benchmark" in payload
    assert "delta" in payload
    assert "ranking_delta_top10_actives" in payload
    assert len(payload["ensemble_rows"]) > 0

    # Vina was called for 2 conformers × 4 compounds
    assert vina_call_count[0] == 8

    # Markdown written
    assert out_md.exists()
    md = out_md.read_text()
    assert "Phase 10" in md
    assert "BEDROC" in md
