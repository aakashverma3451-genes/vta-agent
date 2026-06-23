"""ProteinTTT refinement node tests — hermetic (no GPU, no proteinttt package).

Covers the refinement contract: only low-confidence ESMFold folds are touched,
experimental / AlphaFold / Boltz-2 structures are left alone, and absence or
failure degrades to the unchanged ESMFold fold without crashing.
"""
from __future__ import annotations

import vta.nodes.proteinttt as P
from vta.state import new_state


def _state(structures, proteins=None):
    st = new_state("x.fasta", "ttt-test")
    st["structures"] = structures
    st["extracted_proteins"] = proteins or {}
    return st


def test_refines_low_confidence_esmfold(monkeypatch, tmp_path):
    monkeypatch.setattr(P, "_proteinttt_available", lambda: True)
    monkeypatch.setattr(P, "_refine_structure",
                        lambda seq, init, out: (str(tmp_path / "ref.pdb"), 78.5))

    st = _state(
        {"PB1": {"pdb_path": "/x/PB1.pdb", "mean_plddt": 58.0,
                 "method": "esmfold", "source": "api.esmatlas.com"}},
        {"PB1": {"sequence": "M" * 120}},
    )
    out = P.proteinttt_node(st)
    rec = out["structures"]["PB1"]
    assert rec["method"] == "esmfold+ttt"
    assert rec["mean_plddt"] == 78.5
    assert rec["refined_from_plddt"] == 58.0
    assert rec["pdb_path"] == str(tmp_path / "ref.pdb")
    assert out["versions"]["proteinttt"].startswith("ProteinTTT")
    assert any("58.0 → 78.5" in l for l in out["audit_trail"])


def test_skips_gracefully_when_package_absent(monkeypatch):
    monkeypatch.setattr(P, "_proteinttt_available", lambda: False)

    st = _state({"PB1": {"pdb_path": "/x/PB1.pdb", "mean_plddt": 58.0,
                         "method": "esmfold", "source": "api"}})
    out = P.proteinttt_node(st)
    # Structure untouched; honest skip recorded only because there WAS a candidate.
    assert out["structures"]["PB1"]["method"] == "esmfold"
    assert out["structures"]["PB1"]["mean_plddt"] == 58.0
    assert out["versions"]["proteinttt"] == "skipped"
    assert any("[skip] ProteinTTT" in l for l in out["audit_trail"])


def test_high_confidence_esmfold_is_left_alone(monkeypatch):
    monkeypatch.setattr(P, "_proteinttt_available", lambda: True)
    monkeypatch.setattr(P, "_refine_structure",
                        lambda *a: (_ for _ in ()).throw(AssertionError("must not refine")))

    st = _state({"NS1": {"pdb_path": "/x/NS1.pdb", "mean_plddt": 88.0,
                        "method": "esmfold", "source": "api"}})
    out = P.proteinttt_node(st)
    assert out["structures"]["NS1"]["method"] == "esmfold"     # untouched
    assert "no low-confidence" in out["audit_trail"][-1]
    assert out["versions"]["proteinttt"].startswith("no-op")


def test_experimental_and_alphafold_never_refined(monkeypatch):
    monkeypatch.setattr(P, "_proteinttt_available", lambda: True)
    monkeypatch.setattr(P, "_refine_structure",
                        lambda *a: (_ for _ in ()).throw(AssertionError("must not refine")))

    st = _state({
        # experimental = trusted ground truth, pLDDT None → never a candidate
        "PA":  {"pdb_path": "/x/PA.pdb", "mean_plddt": None,
                "method": "experimental", "source": "8PSO:A"},
        # AlphaFold-DB model, even if low pLDDT, is not an ESMFold fold → skip TTT
        "PB2": {"pdb_path": "/x/PB2.pdb", "mean_plddt": 55.0,
                "method": "alphafold", "source": "AlphaFold DB:P03428"},
    })
    out = P.proteinttt_node(st)
    assert out["structures"]["PA"]["method"] == "experimental"
    assert out["structures"]["PB2"]["method"] == "alphafold"
    assert out["versions"]["proteinttt"].startswith("no-op")


def test_refinement_failure_keeps_esmfold_fold(monkeypatch):
    monkeypatch.setattr(P, "_proteinttt_available", lambda: True)
    monkeypatch.setattr(P, "_refine_structure", lambda *a: None)   # refinement failed

    st = _state(
        {"PB1": {"pdb_path": "/x/PB1.pdb", "mean_plddt": 60.0,
                 "method": "esmfold", "source": "api"}},
        {"PB1": {"sequence": "M" * 100}},
    )
    out = P.proteinttt_node(st)
    rec = out["structures"]["PB1"]
    assert rec["method"] == "esmfold"          # unchanged
    assert rec["mean_plddt"] == 60.0
    assert "refinement failed" in " ".join(out["audit_trail"])
    assert out["versions"]["proteinttt"].startswith("skipped")


def test_noop_when_no_structures():
    out = P.proteinttt_node(new_state("x", "x"))
    assert out.get("versions", {}).get("proteinttt") is None
