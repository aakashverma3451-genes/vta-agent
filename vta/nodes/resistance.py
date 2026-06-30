"""resistance_node — annotation-only resistance-barrier seam."""
from __future__ import annotations

from vta.state import VTAState


def resistance_node(state: VTAState) -> VTAState:
    leads = state.get("lead_candidates") or []
    mutants = state.get("resistance_mutants") or []
    if not leads:
        return state
    if not mutants:
        for lead in leads:
            lead["resistance_barrier"] = {
                "status": "skipped_no_mutant_panel",
                "score": None,
                "ranking_active": False,
            }
        state["versions"]["resistance"] = "skipped"
        state["audit_trail"].append(
            f"[skip] Resistance: no mutant panel supplied for {len(leads)} leads")
        return state

    for lead in leads:
        rows = [m for m in mutants if m.get("ligand") == lead.get("ligand")]
        if not rows:
            score = None
        else:
            penalties = [float(r.get("fitness_cost") or 0.0) - abs(float(r.get("delta_dG") or 0.0))
                         for r in rows]
            score = round(sum(penalties) / len(penalties), 3)
        lead["resistance_barrier"] = {
            "status": "annotated",
            "score": score,
            "mutants_evaluated": len(rows),
            "ranking_active": False,
        }
    state["versions"]["resistance"] = "mutant-panel seam"
    state["audit_trail"].append(
        f"Resistance: annotated {len(leads)} leads from supplied mutant panel")
    return state
