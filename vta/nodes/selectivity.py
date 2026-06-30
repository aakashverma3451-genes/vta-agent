"""selectivity_node — annotation-only counter-target selectivity seam."""
from __future__ import annotations

from vta.state import VTAState


def selectivity_node(state: VTAState) -> VTAState:
    leads = state.get("lead_candidates") or []
    panel = state.get("counter_target_panel") or []
    if not leads:
        return state
    if not panel:
        for lead in leads:
            lead["selectivity"] = {
                "status": "skipped_no_counter_target_panel",
                "score": None,
                "ranking_active": False,
            }
        state["versions"]["selectivity"] = "skipped"
        state["audit_trail"].append(
            f"[skip] Selectivity: no counter-target panel supplied for {len(leads)} leads")
        return state

    for lead in leads:
        primary = abs(float(lead.get("dG") or 0))
        counter_scores = [abs(float(t.get("dG") or 0)) for t in panel if t.get("ligand") == lead.get("ligand")]
        margin = round(primary - max(counter_scores), 3) if counter_scores else None
        lead["selectivity"] = {
            "status": "annotated",
            "score": margin,
            "panel_size": len(counter_scores),
            "ranking_active": False,
        }
    state["versions"]["selectivity"] = "counter-target panel seam"
    state["audit_trail"].append(
        f"Selectivity: annotated {len(leads)} leads from supplied counter-target panel")
    return state
