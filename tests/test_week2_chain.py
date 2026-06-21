"""Week-2 acceptance test: the full Phase-1 chain produces ranked leads.

Proceed path runs structure(stub) → pockets(mock) → dock(mock) → rank(REAL) and
must yield a descending top-20 lead list. Like the spine test this mocks
`classify_genome`, so it is fast and offline; the biology mocks are deterministic
(seeded by run_id), which lets us assert reproducibility.

Run:  PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python -m pytest tests/ -q
"""
from __future__ import annotations

import pytest

pytest.importorskip("langgraph", reason="requires langgraph")

import taxonagent  # noqa: E402
from vta.state import new_state  # noqa: E402


def _contract(confidence: float) -> dict:
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


def test_full_chain_produces_ranked_leads(monkeypatch, mock_structure_net):
    monkeypatch.setattr(taxonagent, "classify_genome", lambda f: _contract(98.0))
    from vta.graph import build_app

    final = build_app().invoke(new_state("seg1.fasta", "week2-run"))

    # PB1 + PA fold experimentally; PB2 (450aa) is refused → pockets on 2 proteins
    assert set(final["structures"]) == {"PB1", "PB2", "PA"}
    assert final["structures"]["PB2"]["method"] == "refused_too_long"
    assert set(final["pockets"]) == {"PB1", "PA"}             # PB2 skipped, no structure

    leads = final["lead_candidates"]
    assert len(leads) == 20                                   # top-20 produced
    scores = [r["score"] for r in leads]
    assert scores == sorted(scores, reverse=True)             # descending
    # records carry the real Vina-shaped fields the ranker scored on
    top = leads[0]
    assert {"protein", "pocket", "ligand", "dG", "rmsd", "le", "score"} <= set(top)

    trail = "\n".join(final["audit_trail"])
    for marker in ("TaxonAgent", "Router", "Pockets", "Docking", "Ranking"):
        assert marker in trail, f"missing audit marker: {marker}"


def test_chain_is_deterministic_for_same_run_id(monkeypatch, mock_structure_net):
    monkeypatch.setattr(taxonagent, "classify_genome", lambda f: _contract(98.0))
    from vta.graph import build_app

    app = build_app()
    a = app.invoke(new_state("seg1.fasta", "same-id"))
    b = app.invoke(new_state("seg1.fasta", "same-id"))
    assert [r["ligand"] for r in a["lead_candidates"]] == \
           [r["ligand"] for r in b["lead_candidates"]]
    assert [r["score"] for r in a["lead_candidates"]] == \
           [r["score"] for r in b["lead_candidates"]]


def test_defer_path_skips_the_whole_chain(monkeypatch):
    monkeypatch.setattr(taxonagent, "classify_genome", lambda f: _contract(40.0))
    from vta.graph import build_app

    final = build_app().invoke(new_state("seg1.fasta", "defer-run"))
    assert final["route"] == "defer"
    assert final["lead_candidates"] == []
    assert final.get("docking_results") is None       # dock node never ran
    assert final.get("pockets") is None               # pocket node never ran
