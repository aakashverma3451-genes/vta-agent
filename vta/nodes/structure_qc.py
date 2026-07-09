"""structure_qc_node — annotate structure quality before pocket detection."""
from __future__ import annotations

from vta.state import VTAState
from vta.structure.qc import qc_structure


def structure_qc_node(state: VTAState) -> VTAState:
    structures = state.get("structures") or {}
    flagged = 0
    for record in structures.values():
        qc = qc_structure(record)
        record["structure_qc"] = qc
        if qc["status"] != "pass":
            flagged += 1
    state["versions"]["structure_qc"] = "pLDDT/resolution/Rfree QC contract"
    state["audit_trail"].append(
        f"Structure QC: annotated {len(structures)} structure(s), {flagged} flagged")
    return state
