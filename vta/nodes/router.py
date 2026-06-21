"""Confidence routing — the mechanism that lets the system say "I don't know".

Plan Task 1.3. This is the single most important reliability feature: without it a
pipeline optimised to produce an answer emits confident output on genuinely
ambiguous input, because nothing surfaces "evidence is insufficient."

IMPORTANT LangGraph note: state mutations made inside a *conditional-edge*
function are discarded — LangGraph uses that function's return value only to pick
the branch. So the decision-and-logging lives in a real NODE (`route_by_confidence`,
whose writes persist), and a separate PURE edge function (`route_edge`) just reads
the recorded decision to choose the next node.

Honesty note (reconciliation mismatch #3): TaxonAgent's `confidence_pct` is an
amino-acid *identity band* (best-hit % identity), not a calibrated probability.
The thresholds come from the Integration Strategy doc, but we read
`classification_basis` so the audit trail states what the number actually means
rather than dressing an identity score up as calibrated confidence.
"""
from __future__ import annotations

from vta.state import VTAState

# Routing thresholds (Integration Strategy doc). Applied to confidence_pct.
PROCEED_THRESHOLD = 95.0     # >= this: proceed cleanly
FLAG_THRESHOLD = 85.0        # [FLAG, PROCEED): proceed but flag for review
# below FLAG_THRESHOLD: defer to a human expert rather than guess

# Node names the edge can route to (constants so the graph wiring and the edge
# function can't disagree on a string).
STRUCTURE_NODE = "structure_node"
DEFER_NODE = "defer_node"


def _basis_phrase(state: VTAState) -> str:
    """Describe the confidence number truthfully for the audit log."""
    basis = (state.get("classification_basis") or "").lower()
    if "identity" in basis:
        return "aa identity"
    if "calibrat" in basis:
        return "calibrated confidence"
    return "confidence"


def route_by_confidence(state: VTAState) -> VTAState:
    """NODE: decide proceed/flag/defer, record it in state + audit trail.

    Writing `route` here (not in the edge) is what makes the decision persist.
    """
    c = state["classification_confidence"]
    unit = _basis_phrase(state)

    if c >= PROCEED_THRESHOLD:
        state["route"] = "proceed"
        state["audit_trail"].append(
            f"Router: {c}% {unit} >= {PROCEED_THRESHOLD:g} -> PROCEED")
    elif c >= FLAG_THRESHOLD:
        state["route"] = "flag"
        state["audit_trail"].append(
            f"Router: {c}% {unit} in [{FLAG_THRESHOLD:g},{PROCEED_THRESHOLD:g}) "
            f"-> PROCEED w/ FLAG")
    else:
        state["route"] = "defer"
        state["audit_trail"].append(
            f"Router: {c}% {unit} < {FLAG_THRESHOLD:g} -> DEFER to expert")
    return state


def route_edge(state: VTAState) -> str:
    """PURE conditional edge: read the recorded decision, return the next node."""
    return DEFER_NODE if state["route"] == "defer" else STRUCTURE_NODE
