"""Week-1 acceptance test: the spine runs end-to-end.

Proves the *wiring* (classify → route → structure|defer → END + audit trail),
not TaxonAgent itself. So we monkeypatch `taxonagent.classify_genome` with a
canned §4.1.1 contract — the spine test stays fast and offline (no DIAMOND DB,
no network), and TaxonAgent has its own suite for the classification itself.

Run:  PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python -m pytest tests/ -q
"""
from __future__ import annotations

import pytest

pytest.importorskip("langgraph", reason="Task 1.4 requires langgraph")

import taxonagent  # noqa: E402
from vta.state import new_state  # noqa: E402


def _contract(confidence: float) -> dict:
    """A minimal but schema-valid §4.1.1 contract for a TiLV-like input."""
    return {
        "strain_id": "amnoonviridae_1",
        "genus": "Tilapinevirus",
        "species": "Tilapia lake virus",
        "confidence_pct": confidence,
        "confidence_basis": "best-hit amino-acid identity (DIAMOND blastp) — "
                            "an identity band, not a calibrated probability",
        "ictv_rule_match": [],
        "extracted_proteins": {
            "PB1": {"sequence": "M" * 560, "length_aa": 560, "plddt": None},
            "PB2": {"sequence": "K" * 450, "length_aa": 450, "plddt": None},
            "PA":  {"sequence": "L" * 340, "length_aa": 340, "plddt": None},
        },
        "reference_alignment": {},
        "audit_trail": [],
        "timestamp": "2026-06-22T00:00:00Z",
        "taxonagent_version": taxonagent.__version__,
        "kg_version": "ICTV-VMR-MSL40",
    }


@pytest.fixture
def app(monkeypatch):
    # Patch the real entry point with a canned contract → no DIAMOND, no network.
    from vta.graph import build_app
    return build_app


def test_spine_proceeds_on_high_confidence(monkeypatch, mock_structure_net):
    monkeypatch.setattr(taxonagent, "classify_genome", lambda f: _contract(98.0))
    from vta.graph import build_app

    final = build_app().invoke(new_state("seg1.fasta", "spine-proceed"))

    assert final["classification_confidence"] == 98.0
    assert final["route"] == "proceed"
    assert set(final["structures"]) == {"PB1", "PB2", "PA"}   # structure node ran
    trail = "\n".join(final["audit_trail"])
    assert "TaxonAgent" in trail and "PROCEED" in trail and "Structure[" in trail


def test_spine_defers_on_low_confidence(monkeypatch):
    monkeypatch.setattr(taxonagent, "classify_genome", lambda f: _contract(40.0))
    from vta.graph import build_app

    final = build_app().invoke(new_state("seg1.fasta", "spine-defer"))

    assert final["route"] == "defer"
    assert final["lead_candidates"] == []        # defer node ran, folded nothing
    assert final.get("structures") is None       # structure node did NOT run
    assert "DEFER" in "\n".join(final["audit_trail"])
