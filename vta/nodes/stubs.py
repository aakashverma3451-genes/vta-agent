"""Placeholder nodes that let the graph run end-to-end before the real biology.

Plan Task 1.4. The structure stub lets the Week-1 spine compile and execute —
genome → classify → route → (structure | defer) → END — proving the wiring before
ESMFold (Task 2.1) replaces `structure_node_stub` with the real folder.

`defer_node` is NOT a stub: it is the real, permanent honest-exit node used when
the router decides confidence is too low to proceed.
"""
from __future__ import annotations

from vta.state import VTAState


def structure_node_stub(state: VTAState) -> VTAState:
    """Pass-through stand-in for the real ESMFold node. Logs intent, fills shape."""
    proteins = state.get("extracted_proteins") or {}
    state["audit_trail"].append(
        f"STUB structure_node: would fold {list(proteins)}"
    )
    state["structures"] = {
        name: {"pdb_path": None, "mean_plddt": None} for name in proteins
    }
    return state


def defer_node(state: VTAState) -> VTAState:
    """Honest exit: confidence too low — return hypotheses to a human, fold nothing."""
    state["audit_trail"].append(
        "DEFER: returning top-3 genus hypotheses to expert (no structures folded)"
    )
    state["lead_candidates"] = []
    return state
