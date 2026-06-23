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


# ── AlphaFold DB cascade (>400 aa with a UniProt accession) ───────────────────
def test_alphafold_used_for_large_protein_with_accession(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    # >400 aa, not in EXPERIMENTAL_PDB, but carries a UniProt accession -> AlphaFold DB
    # (no length cap) instead of refusing. Boltz-2 is available but must NOT be used.
    monkeypatch.setattr(S, "fetch_alphafold", lambda acc: fake_esmfold_pdb(plddt=91.0))
    monkeypatch.setattr(S, "_boltz2_bin", lambda: "boltz")
    st = _state({"PB2": {"sequence": "K" * 450, "uniprot": "P03428"}})
    out = S.structure_node(st)
    rec = out["structures"]["PB2"]
    assert rec["method"] == "alphafold"
    assert rec["mean_plddt"] == 91.0
    assert rec["source"] == "AlphaFold DB:P03428"
    assert any("AlphaFold DB P03428" in line for line in out["audit_trail"])


def test_uniprot_query_resolves_accession_then_alphafold(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    # >400 aa, no explicit accession, but an opt-in free-text query -> resolve via
    # UniProt, then fold via AlphaFold DB. Resolver + AlphaFold are both seams.
    monkeypatch.setattr(S, "resolve_uniprot", lambda q: "P03431")
    monkeypatch.setattr(S, "fetch_alphafold", lambda acc: fake_esmfold_pdb(plddt=88.0))
    st = _state({"NSP12": {"sequence": "M" * 600,    # novel name, not in EXPERIMENTAL_PDB
                           "uniprot_query": "NSP12 SARS-CoV-2 RdRp"}})
    out = S.structure_node(st)
    rec = out["structures"]["NSP12"]
    assert rec["method"] == "alphafold"
    assert rec["source"] == "AlphaFold DB:P03431"
    assert any("resolved UniProt P03431" in line for line in out["audit_trail"])


def test_uniprot_query_unresolved_falls_through_to_refuse(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    # Resolver finds nothing and no Boltz-2 -> honest refusal, no crash.
    monkeypatch.setattr(S, "resolve_uniprot", lambda q: None)
    monkeypatch.setattr(S, "_boltz2_bin", lambda: None)
    st = _state({"NSP12": {"sequence": "M" * 600, "uniprot_query": "gibberish"}})
    out = S.structure_node(st)
    assert out["structures"]["NSP12"]["method"] == "refused_too_long"


def test_alphafold_404_degrades_gracefully(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    def not_found(acc):
        raise RuntimeError("HTTP 404")
    monkeypatch.setattr(S, "fetch_alphafold", not_found)
    st = _state({"PB2": {"sequence": "K" * 450, "accession": "Q99999"}})
    out = S.structure_node(st)                         # must NOT raise
    rec = out["structures"]["PB2"]
    assert rec["method"] == "alphafold_failed"
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
