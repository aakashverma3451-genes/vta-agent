"""Supplemental report sections kept separate from the core renderer."""
from __future__ import annotations

from typing import Any

from vta.provenance import run_provenance
from vta.state import VTAState


def scientific_status(state: VTAState, esc) -> str:
    leads = state.get("lead_candidates") or []
    if not leads:
        return ""
    rows = []
    for lead in leads[:10]:
        sel = lead.get("selectivity") or {}
        res = lead.get("resistance_barrier") or {}
        rows.append(
            f"<tr><td>{esc(lead.get('ligand'))}</td>"
            f"<td>{esc(sel.get('score'))}</td><td>{esc(sel.get('status'))}</td>"
            f"<td>{esc(res.get('score'))}</td><td>{esc(res.get('status'))}</td></tr>"
        )
    registry = state.get("prediction_registry") or {}
    head = "<tr><th>Ligand</th><th>Selectivity</th><th>Status</th><th>Resistance</th><th>Status</th></tr>"
    note = (
        "<p class='note'><b>Scientific status</b>: outputs are ranked hypotheses until "
        "prospectively locked and experimentally confirmed. Research use only; do not "
        "use this system for pathogen enhancement.</p>"
    )
    reg = f"<p class='note'>Prediction registry: {esc(registry.get('status') or 'not locked')}</p>"
    return f"<h2>Scientific validation status</h2><table>{head}{''.join(rows)}</table>{reg}{note}"


def footer_provenance(state: VTAState, esc) -> str:
    return f"<h2>Run provenance</h2><pre>{esc(run_provenance(state))}</pre>"
