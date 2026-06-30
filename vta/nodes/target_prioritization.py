"""target_prioritization_node — antiviral target ranking seam."""
from __future__ import annotations

from vta.state import VTAState

WEIGHTS = {
    "druggability": 0.25,
    "resistance_liability": 0.20,
    "essentiality": 0.20,
    "host_selectivity": 0.15,
    "structural_confidence": 0.10,
    "literature_evidence": 0.10,
}


def _score(record: dict) -> float:
    total = 0.0
    for key, weight in WEIGHTS.items():
        value = float(record.get(key, 0.5))
        if key == "resistance_liability":
            value = 1.0 - value
        total += weight * max(0.0, min(1.0, value))
    return round(total, 4)


def target_prioritization_node(state: VTAState) -> VTAState:
    targets = state.get("candidate_targets")
    if targets is None:
        proteins = state.get("extracted_proteins") or {}
        targets = [{"target": name, "structural_confidence": (rec.get("plddt") or 50) / 100}
                   for name, rec in proteins.items()]
    prioritized = []
    for target in targets:
        row = dict(target)
        row["target_priority_score"] = _score(row)
        row["conservation_role"] = "target_prioritization"
        prioritized.append(row)
    prioritized.sort(key=lambda r: r["target_priority_score"], reverse=True)
    state["target_prioritization"] = prioritized
    state["versions"]["target_prioritization"] = "weighted antiviral target criteria v1"
    state["audit_trail"].append(
        f"Target prioritization: scored {len(prioritized)} target(s); compound ranking unchanged")
    return state
