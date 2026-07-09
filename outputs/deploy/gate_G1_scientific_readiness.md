# Gate G1 — Scientific Readiness (Deployment plan D0)

**Date:** 2026-07-03
**Decision:** **PASS** — the D0 scientific-readiness gate is met. Service work (D1+) may begin.
**Owner:** session-a (opus)
**Scope:** D0.1–D0.5 of `VTA-Agent_Deployment_Implementation_Plan.md`. This gate blocks all
service/MLOps work; nothing in D1–D10 was started before this decision.

> Acceptance framing (from the plan): *"the report is honest by construction and survives peer
> review," not "tests green."* Every ranked claim a deployed user sees must be bounded by a
> benchmark grade + CI, docking must be shown to beat-or-not-beat a trivial baseline by a
> paired test, pose reliability must be stated, ranking weights must have a validation basis,
> and the hypotheses-not-efficacy framing must be present and non-removable.

---

## How each D0 criterion is met

### D0.1 — Trivial-baseline floor, embedded in the report — **MET**
- Baseline engine: `vta/eval/baselines.py` (ECFP4 leave-one-out 2-D similarity + random null),
  committed in Phase 11 (WI-1).
- **Now wired into the output**: the honesty envelope (`vta/report_envelope.py`) reads the
  committed paired verdict from `outputs/phase11/baselines_and_gate.json` and renders, per
  target, *"Docking vs trivial baseline: beats / does NOT beat 2D-sim"*. On the powered Mpro
  benchmark the rendered verdict is the honest negative: **Vina does NOT beat 2D-similarity**
  (`structure_based_enrichment_demonstrated = false`).

### D0.2 — Redocking / pose sanity, gated — **MET**
- Redocking check: `outputs/phase11/redock_validation.json` (Phase 11 WI-2).
- **Now surfaced per target** in the envelope: Mpro 7L11/XF1 → **pose-reliable (1.654 Å < 2.0)**;
  HCV NS5B 2XI3/GTP → **pose reliability UNKNOWN (RMSD uncomputed)** — labelled, not hidden.

### D0.3 — Paired significance, not CI-overlap — **MET**
- `vta/eval/significance.py` (`paired_bootstrap_delta`, `holm_bonferroni`); all CI-overlap
  significance logic was purged from the codebase in Phase 11 WI-3 (incl. the Phase-10 scripts).
- Every promote/do-not-promote decision the envelope cites (LE gate, baseline verdict) is a
  paired-bootstrap result, Holm-corrected across the metric family.

### D0.4 — Ligand-efficiency demotion — **MET**
- `vta/nodes/rank.py`: `WEIGHTS = {"le": 0.0, "dG": 1.0, "conservation": 0.0}` — ΔG-primary.
  LE and conservation are reported annotations (weight 0). This was decided by the WI-6 paired
  gate: the LE-led composite did not beat ΔG-only on the powered Mpro benchmark
  (paired Δ BEDROC 95% CI [−0.58, +0.06]; `live_rank_weights` recorded in
  `outputs/phase11/baselines_and_gate.json`).
- **Stale-honesty defects fixed this pass** (the report previously *misdescribed its own
  ranking*): the leads note that said *"scored by ligand efficiency, not raw ΔG"* and the two
  conservation notes calling conservation *"a 10% term in the score"* were rewritten to the
  in-force ΔG-primary / weight-0-annotation reality. Stale provenance label
  `"LE-led composite v1"` → `"ΔG-primary (LE demoted, Phase 11 WI-6 gate)"`. No weight changed.

### D0.5 — Report-level honesty contract — **MET (new work this pass)**
- New `vta/report_envelope.py` : `build_envelope(state)` assembles the non-removable envelope
  from committed artifacts only (no inline-asserted numbers):
  - **disclaimer** — "ranked, uncertainty-bearing HYPOTHESES … NOT clinical, efficacy, or
    safety claims";
  - **benchmark pin** — artifact + `frozen_at` + `content_hash` (`fdac0664fe0cea8f`) +
    gate decision, from `outputs/phase9/locked_benchmark.json` (D3 benchmark-pinning link);
  - **per-target grade + metric 95% CI**, trivial-baseline verdict, pose reliability;
  - **live ranking weights** imported from `vta.nodes.rank` (can never drift from what ranked);
  - **scoring caveats** — the Phase-9D nucleotide/metal caveat + the frozen benchmark's own
    modest-signal / underpowered caveats;
  - **out-of-validated-domain flag** for any run target no frozen benchmark covers.
- **Structural, non-removable:** `render_report` calls `_envelope()` unconditionally as the
  first card, and `report_node` stores the same object on `state["honesty_envelope"]` (so the
  future API response and run manifest read one source). Tests assert the disclaimer +
  benchmark context appear for a full run, a deferred run, and an empty run — there is no code
  path that emits a ranking without the envelope.

---

## Evidence
- Code: `vta/report_envelope.py`, `vta/report.py` (envelope card + corrected notes),
  `vta/nodes/report.py`, `vta/state.py` (`honesty_envelope`), `vta/nodes/rank.py` (label).
- Tests: `tests/test_report_envelope.py` (hermetic; builder behaviour + non-removable
  guarantee) + updated `tests/test_report.py`. **Full suite: 240 passed.**
- Sourced artifacts (read, not edited): `outputs/phase9/locked_benchmark.json`,
  `outputs/phase11/baselines_and_gate.json`, `outputs/phase11/redock_validation.json`,
  `outputs/phase11/mpro_30to1_benchmark.json`.

## Honest limitations recorded at this gate (not blockers, but carried forward)
- The powered Mpro docking signal is **modest and does not beat a trivial 2-D baseline**; the
  envelope states this on every Mpro report rather than hiding it. Deployment does not improve
  the signal — it prevents the signal from being over-read.
- HCV NS5B pose reliability is **unknown** (triphosphate RMSD tooling gap); the envelope says so.
- The out-of-validated-domain flag is keyword-based (Mpro / NS5B-HCV / PB1-TiLV). Any other
  target is correctly flagged exploratory, but the map must be extended as targets are added.

## What G1 does NOT authorize
- No claim that VTA-Agent's docking has demonstrated structure-based skill (it has not).
- No ranking-weight or rescorer change — those remain behind the D0.3 paired-test gate.
- No service, persistence, container, or governance work is validated by this gate; those are
  G2/G3 and are not yet started.
