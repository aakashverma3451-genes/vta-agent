"""Unit tests for the real structure node (Task 2.1).

Covers the four decisions: experimental-first + chain extraction, ESMFold within
the ceiling, honest refusal above it, and graceful failure on a network error.
Pure helpers (extract_chain, mean_plddt) are tested directly. Network is mocked.
"""
from __future__ import annotations

import pytest

import vta.nodes.structure as S
from tests.conftest import fake_complex_pdb, fake_esmfold_pdb
from vta.state import new_state


def _state(proteins: dict) -> dict:
    st = new_state("x.fasta", "struct-test")
    st["extracted_proteins"] = proteins
    return st


# ── pure helpers ─────────────────────────────────────────────────────────────
def test_extract_chain_keeps_only_that_chain():
    pdb = fake_complex_pdb(chains=("A", "B"), n_res=10)
    a = S.extract_chain(pdb, "A")
    assert a.count(" CA ") == 10
    assert all(line[21] == "A" for line in a.splitlines() if line.startswith("ATOM"))
    assert "ATOM" in a and a.rstrip().endswith("END")


def test_mean_plddt_reads_bfactor_column():
    assert S.mean_plddt(fake_esmfold_pdb(plddt=88.0, n_res=30)) == 88.0


# ── node decisions ───────────────────────────────────────────────────────────
def test_experimental_first_extracts_chain_and_keeps_plddt_none(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(S, "fetch_rcsb_pdb", lambda pid: fake_complex_pdb(("A", "B")))

    st = _state({"PB1": {"sequence": "M" * 500}, "PA": {"sequence": "L" * 316}})
    out = S.structure_node(st)

    for name in ("PB1", "PA"):
        rec = out["structures"][name]
        assert rec["method"] == "experimental"
        assert rec["mean_plddt"] is None         # trusted ground truth, NOT low
        assert rec["pdb_path"] and rec["source"].startswith("8PSO:")


def test_pb2_over_ceiling_is_refused_not_crashed(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    # PB2 has no experimental chain and is >400 aa -> must refuse, not fold/crash.
    st = _state({"PB2": {"sequence": "K" * 450}})
    out = S.structure_node(st)
    rec = out["structures"]["PB2"]
    assert rec["method"] == "refused_too_long"
    assert rec["pdb_path"] is None
    assert any("REFUSED" in line for line in out["audit_trail"])


def test_esmfold_used_for_small_novel_protein(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(S, "fold_esmfold", lambda seq: fake_esmfold_pdb(plddt=82.0))
    st = _state({"NS1": {"sequence": "Q" * 120}})     # not in EXPERIMENTAL_PDB, <400
    out = S.structure_node(st)
    rec = out["structures"]["NS1"]
    assert rec["method"] == "esmfold"
    assert rec["mean_plddt"] == 82.0


def test_experimental_fetch_failure_degrades_gracefully(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    def boom(pdb_id):
        raise RuntimeError("network down")
    monkeypatch.setattr(S, "fetch_rcsb_pdb", boom)

    st = _state({"PB1": {"sequence": "M" * 500}})
    out = S.structure_node(st)                         # must NOT raise
    rec = out["structures"]["PB1"]
    assert rec["method"] == "experimental_failed"
    assert rec["pdb_path"] is None


# ── Boltz-2 cascade ──────────────────────────────────────────────────────────
def test_boltz2_used_for_large_protein_when_available(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(S, "_boltz2_bin", lambda: "boltz")
    monkeypatch.setattr(S, "fold_boltz2",
                        lambda seq, name, out_dir: (str(tmp_path / f"{name}.pdb"), 74.5))

    st = _state({"PB2": {"sequence": "K" * 450}})  # >400 aa, no experimental chain
    out = S.structure_node(st)
    rec = out["structures"]["PB2"]
    assert rec["method"] == "boltz2"
    assert rec["mean_plddt"] == 74.5
    assert rec["pdb_path"] is not None
    assert any("Boltz-2" in l for l in out["audit_trail"])
    assert "Boltz-2" in out["versions"]["esmfold"]


def test_boltz2_failure_still_records_no_structure(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(S, "_boltz2_bin", lambda: "boltz")
    monkeypatch.setattr(S, "fold_boltz2", lambda *a, **kw: None)  # returns None

    st = _state({"PB2": {"sequence": "K" * 450}})
    out = S.structure_node(st)
    assert out["structures"]["PB2"]["method"] == "boltz2_failed"
    assert out["structures"]["PB2"]["pdb_path"] is None


def test_boltz2_not_used_when_absent(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(S, "_boltz2_bin", lambda: None)

    st = _state({"PB2": {"sequence": "K" * 450}})
    out = S.structure_node(st)
    # Without Boltz-2 it should fall back to refused_too_long
    assert out["structures"]["PB2"]["method"] == "refused_too_long"


def test_esmfold_still_used_below_ceiling_even_when_boltz2_present(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(S, "_boltz2_bin", lambda: "boltz")
    monkeypatch.setattr(S, "fold_esmfold", lambda seq: fake_esmfold_pdb(plddt=80.0))

    st = _state({"NS1": {"sequence": "Q" * 120}})   # <400 aa → ESMFold wins
    out = S.structure_node(st)
    assert out["structures"]["NS1"]["method"] == "esmfold"
