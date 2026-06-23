"""md_select_node — choose the top-N ranked leads for MD validation.

The first node of the opt-in MD phase (Phase 4). Takes the ADMET-annotated
lead_candidates, selects the top MD_TOP_N by composite score, and writes
them into `md_candidates` so md_simulate has a bounded, labelled work list.

Runs on CPU, no external deps.
"""
from __future__ import annotations

from vta.state import VTAState

MD_TOP_N = 5


def md_select_node(state: VTAState) -> VTAState:
    leads = state.get("lead_candidates") or []
    if not leads:
        state["md_candidates"] = []
        state["audit_trail"].append("MD-select: no leads — MD phase skipped")
        return state
    selected = leads[:MD_TOP_N]
    state["md_candidates"] = selected
    names = [f"{s['ligand']} ({s['score']})" for s in selected]
    state["audit_trail"].append(
        f"MD-select: {len(selected)} leads queued for MD validation: "
        + ", ".join(names)
    )
    return state
