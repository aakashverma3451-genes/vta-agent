"""consensus_node — annotation-only agreement across orthogonal scoring seams."""
from __future__ import annotations

from vta.state import VTAState


def _norm(values: list[float], value: float, invert: bool = False) -> float:
    if not values:
        return 0.0
    lo, hi = min(values), max(values)
    if hi == lo:
        return 0.5
    scaled = (value - lo) / (hi - lo)
    return round(1.0 - scaled if invert else scaled, 4)


def consensus_node(state: VTAState) -> VTAState:
    rows = state.get("docking_results") or []
    if not rows:
        return state
    dgs = [float(r["dG"]) for r in rows if r.get("dG") is not None]
    cnn = [float(r["cnn_affinity"]) for r in rows if r.get("cnn_affinity") is not None]
    boltz = [float(r["boltzina_score"]) for r in rows if r.get("boltzina_score") is not None]
    annotated = 0
    for r in rows:
        signals = {"vina": _norm(dgs, float(r["dG"]), invert=True)} if r.get("dG") is not None else {}
        if r.get("cnn_affinity") is not None:
            signals["gnina"] = _norm(cnn, float(r["cnn_affinity"]))
        if r.get("boltzina_score") is not None:
            signals["boltzina"] = _norm(boltz, float(r["boltzina_score"]))
        if signals:
            r["consensus"] = {
                "score": round(sum(signals.values()) / len(signals), 4),
                "signals": signals,
                "ranking_active": False,
            }
            annotated += 1
    state["versions"]["consensus"] = "annotation-only vina/gnina/boltzina agreement"
    state["audit_trail"].append(
        f"Consensus: annotated {annotated} docking records (annotation-only, ranking unchanged)")
    return state
