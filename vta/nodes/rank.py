"""rank_node (REAL LOGIC) — turn docking records into ranked lead candidates.

Plan Task 2.3 (the ranking half). This is NOT mocked: it encodes the scientific
judgment (which criteria matter and how much), so it is built for real now and is
unchanged when real Vina replaces the mock in Phase 2.

Ranking is by **AutoDock Vina affinity (ΔG)** — the primary and only ranking term.
Ligand efficiency (LE) and conservation are computed and REPORTED per compound (on every
record) but are NOT in the ranking.

Why this changed — a GATED decision, not arbitrary (Phase 11 WI-6):
Earlier versions led with LE (0.55 LE / 0.35 ΔG / 0.10 conservation), justified by a
9-ligand control-recovery check that appeared to fail under ΔG-led weights. Phase 8 showed
that check was INFLATED — control recovery collapsed 3/4 → 1/4 once property-matched decoys
were used. The Phase 11 gate (`scripts/phase11_baselines_gate.py`) then tested the LE-led
composite against ΔG-only on the powered Mpro benchmark with a PAIRED bootstrap (not CI
overlap): the composite did NOT beat ΔG-only on BEDROC (median Δ = −0.27, 95% CI
[−0.58, +0.06]). Per the WI-6 rule, LE is therefore demoted to a reported secondary
annotation and ΔG becomes the primary ranking term. This is consistent with Kenny 2019
(J. Cheminform. 11:8): LE is size- and unit-dependent and has no benchmark basis for leading
a ranking.

Normalisation is min/max once per criterion, applied inline; `le` and `conservation` remain
on each record as reported annotations (weight 0 in the score).
"""
from __future__ import annotations

from vta.provenance import score_provenance
from vta.state import VTAState

# ΔG-primary after the Phase 11 WI-6 gate: on the powered Mpro benchmark the LE-led composite
# did not beat ΔG-only (paired-bootstrap Δ BEDROC 95% CI [−0.58, +0.06]). LE and conservation
# are reported per-compound annotations (weight 0), not ranking terms. Change is gated, with
# the supporting Δ-CI logged in outputs/phase11/ (do not revert without a new gate result).
WEIGHTS = {"le": 0.0, "dG": 1.0, "conservation": 0.0}
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
    le_good = [-r["le"] for r in rows]        # more negative le -> higher goodness
    dG_lo, dG_hi = _minmax(dG_good)
    le_lo, le_hi = _minmax(le_good)

    for r in rows:
        r["score"] = round(
            WEIGHTS["le"] * _norm(-r["le"], le_lo, le_hi)
            + WEIGHTS["dG"] * _norm(-r["dG"], dG_lo, dG_hi)
            + WEIGHTS["conservation"] * r["conservation"],
            4,
        )
        r["score_provenance"] = score_provenance(
            "VTA rank_node", "LE-led composite v1", inputs={
                "weights": WEIGHTS,
                "ligand": r.get("ligand"),
                "protein": r.get("protein"),
                "pocket": r.get("pocket"),
                "dG": r.get("dG"),
                "le": r.get("le"),
                "conservation": r.get("conservation"),
            })

    top = sorted(rows, key=lambda r: r["score"], reverse=True)[:TOP_N]
    state["lead_candidates"] = top
    best = top[0]
    state["audit_trail"].append(
        f"Ranking: top lead {best['ligand']} on {best['protein']} "
        f"score {best['score']} (dG {best['dG']}, conservation {best['conservation']})"
    )
    return state
