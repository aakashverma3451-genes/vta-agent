"""md_rerank_node — re-rank leads with MD stability evidence (Phase 4).

Assigns an MD verdict (STABLE / MODERATE / UNSTABLE / not run) to each
md_candidate and produces `md_validated_leads` — the Phase-4 output that
the report renders as its "MD-validated" table.

Score adjustment philosophy (conservative by design):
  - MD should CONFIRM good docking leads, not surface new ones.
  - STABLE poses get a 20% boost — independent vote for the lead.
  - MODERATE poses get a 10% penalty — ligand moves but stays in pocket.
  - UNSTABLE poses get a 70% penalty — they drifted; discard unless the
    docking score was exceptional.
  These multipliers are the same ones used in the MD Implementation Guide
  and in Yamaotsu & Hirono 2016 (the reference validation methodology).

The composite score here uses the docking `score` field set by rank_node
(LE-led). We do NOT re-open the rank_node weighting; the MD term is a
post-hoc overlay that flags confidence, not a new scoring dimension.
"""
from __future__ import annotations

from vta.state import VTAState

_BOOST_STABLE = 1.20
_PENALTY_MODERATE = 0.90
_PENALTY_UNSTABLE = 0.30


def md_rerank_node(state: VTAState) -> VTAState:
    candidates = state.get("md_candidates") or []
    analysis = state.get("md_analysis") or {}

    for lead in candidates:
        lid = lead["ligand"]
        verdict = (analysis.get(lid) or {}).get("md_verdict", "not_run")
        base = lead.get("score", 0.0)

        if verdict == "STABLE":
            lead["md_score"] = round(base * _BOOST_STABLE, 4)
            lead["md_badge"] = "MD-STABLE"
        elif verdict == "MODERATE":
            lead["md_score"] = round(base * _PENALTY_MODERATE, 4)
            lead["md_badge"] = "MD-MODERATE"
        elif verdict == "UNSTABLE":
            lead["md_score"] = round(base * _PENALTY_UNSTABLE, 4)
            lead["md_badge"] = "MD-UNSTABLE"
        else:
            lead["md_score"] = base
            lead["md_badge"] = "MD not run"

        lead["md_details"] = analysis.get(lid) or {}

    ranked = sorted(candidates, key=lambda x: x.get("md_score", 0), reverse=True)
    state["md_validated_leads"] = ranked
    top = ranked[0] if ranked else {}
    state["audit_trail"].append(
        f"MD-rerank: top MD-validated lead = {top.get('ligand', 'none')} "
        f"({top.get('md_badge', '')}, md_score {top.get('md_score')})"
    )
    return state
