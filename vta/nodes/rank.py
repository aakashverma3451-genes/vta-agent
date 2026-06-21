"""rank_node (REAL LOGIC) — turn docking records into ranked lead candidates.

Plan Task 2.3 (the ranking half). This is NOT mocked: it encodes the scientific
judgment (which criteria matter and how much), so it is built for real now and is
unchanged when real Vina replaces the mock in Phase 2.

Composite score per ligand-pocket record, each criterion min-max normalised to
[0,1] where 1 = best, then weighted:

    score = 0.50 * dG'        (binding affinity — more negative is better)
          + 0.25 * rmsd'      (pose stability — lower is better)
          + 0.15 * le'        (ligand efficiency — more negative per heavy atom)
          + 0.10 * conservation (target conservation — higher is better)

Implementation note (fixes reconciliation snag #4): normalisation is computed from
each criterion's min/max ONCE and applied to each row inline — no dict keyed by
float value (which collides silently and is fragile). `conservation` is already on
a 0–1 scale, so it is used directly.
"""
from __future__ import annotations

from vta.state import VTAState

WEIGHTS = {"dG": 0.50, "rmsd": 0.25, "le": 0.15, "conservation": 0.10}
TOP_N = 20


def _minmax(values: list[float]) -> tuple[float, float]:
    return min(values), max(values)


def _norm(val: float, lo: float, hi: float) -> float:
    """Min-max to [0,1]; neutral 1.0 when the criterion is constant (no signal)."""
    return 1.0 if hi == lo else (val - lo) / (hi - lo)


def rank_node(state: VTAState) -> VTAState:
    rows = state.get("docking_results") or []
    if not rows:
        state["lead_candidates"] = []
        state["audit_trail"].append("Ranking: no docking results to rank")
        return state

    # "Goodness" transforms so higher == better for every criterion, then one
    # min/max per criterion (computed once, not per-row dict lookup).
    dG_good = [-r["dG"] for r in rows]        # more negative dG -> higher goodness
    rmsd_good = [-r["rmsd"] for r in rows]    # lower rmsd -> higher goodness
    le_good = [-r["le"] for r in rows]        # more negative le -> higher goodness
    dG_lo, dG_hi = _minmax(dG_good)
    rmsd_lo, rmsd_hi = _minmax(rmsd_good)
    le_lo, le_hi = _minmax(le_good)

    for r in rows:
        r["score"] = round(
            WEIGHTS["dG"] * _norm(-r["dG"], dG_lo, dG_hi)
            + WEIGHTS["rmsd"] * _norm(-r["rmsd"], rmsd_lo, rmsd_hi)
            + WEIGHTS["le"] * _norm(-r["le"], le_lo, le_hi)
            + WEIGHTS["conservation"] * r["conservation"],
            4,
        )

    top = sorted(rows, key=lambda r: r["score"], reverse=True)[:TOP_N]
    state["lead_candidates"] = top
    best = top[0]
    state["audit_trail"].append(
        f"Ranking: top lead {best['ligand']} on {best['protein']} "
        f"score {best['score']} (dG {best['dG']}, conservation {best['conservation']})"
    )
    return state
