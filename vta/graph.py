"""vta.graph — compile the VTA-Agent LangGraph application.

Plan Task 1.4 (Week-1 spine). Wires the proven nodes into a runnable graph:

    classify ──▶ route_by_confidence ──┬─(proceed/flag)─▶ structure_node ─▶ END
                                       └─(defer)────────▶ defer_node     ─▶ END

`structure_node` is the stub here; Task 2.4 swaps in the real ESMFold node and
appends the pockets → dock → rank chain. The router returns the next node's name,
so the conditional-edge mapping below must list every name `route_by_confidence`
can return (STRUCTURE_NODE / DEFER_NODE).
"""
from __future__ import annotations

from langgraph.graph import END, StateGraph

from vta.nodes.classify import classify_node
from vta.nodes.docking import docking_node
from vta.nodes.pockets import pockets_node
from vta.nodes.admet import admet_node
from vta.nodes.rank import rank_node
from vta.nodes.rescore import rescore_node
from vta.nodes.report import report_node
from vta.nodes.router import (
    DEFER_NODE,
    STRUCTURE_NODE,
    route_by_confidence,
    route_edge,
)
from vta.nodes.stubs import defer_node
from vta.nodes.structure import structure_node
from vta.state import VTAState


def build_app():
    """Compile the full graph: classify → route → (chain | defer) → report → END.

    Proceed path:  structure → pockets → dock → rescore → rank → admet ─┐
    Defer path:    defer ──────────────────────────────────────────────┴─▶ report → END

    Both paths end at report_node, so every run emits a self-contained HTML report
    (a deferred run's report shows the classification + defer reason, no leads).
    """
    g = StateGraph(VTAState)

    g.add_node("classify", classify_node)
    g.add_node("router", route_by_confidence)      # NODE: writes `route` (persists)
    g.add_node(STRUCTURE_NODE, structure_node)
    g.add_node("pockets", pockets_node)
    g.add_node("dock", docking_node)
    g.add_node("rescore", rescore_node)    # DL CNN re-score seam (annotation-only)
    g.add_node("rank", rank_node)
    g.add_node("admet", admet_node)
    g.add_node(DEFER_NODE, defer_node)
    g.add_node("report", report_node)

    g.set_entry_point("classify")
    g.add_edge("classify", "router")
    # PURE edge reads the recorded decision; keys must cover what route_edge returns.
    g.add_conditional_edges("router", route_edge, {
        STRUCTURE_NODE: STRUCTURE_NODE,
        DEFER_NODE: DEFER_NODE,
    })
    g.add_edge(STRUCTURE_NODE, "pockets")
    g.add_edge("pockets", "dock")
    g.add_edge("dock", "rescore")      # DL CNN re-score of poses (seam; skips w/o gnina)
    g.add_edge("rescore", "rank")
    g.add_edge("rank", "admet")        # autonomous drug-likeness/safety annotation
    g.add_edge("admet", "report")
    g.add_edge(DEFER_NODE, "report")
    g.add_edge("report", END)

    return g.compile()
