# Phase R — Stage 5 gate (R6 component ablation): PASS

**Date:** 2026-07-07 · **Owner:** session-a (opus)

**Go criterion (from plan):** the ablation shows the workflow layers improve **decision quality
and honesty** even though raw docking does not beat the 2D baseline — i.e. the *architecture*,
not the scorer, is the contribution, with the gain **routing-type-dependent** (HemaGuide's
"no single component sufficient across all case types").

## Result (`scripts/phaseR_ablation.py` → `outputs/phaseR/ablation.{json,md}`)

Held-out targets spanning classes, scored on the emitted CLAIM vs a curated ground-truth
disposition (NOT enrichment):

| Level | decision accuracy | false-confidence rate | OOD-flag recall |
|---|---|---|---|
| L0 raw (rigid Vina + rank) | 0.0 | 1.0 | 0.0 |
| L1 +dossier | 0.0 | 1.0 | 0.0 |
| L2 +triage | 0.6 | 0.4 | 1.0 |
| L3 +verification | 1.0 | 0.0 | 1.0 |
| L4 full agent | 1.0 | 0.0 | 1.0 |

**Which layer fixed which target:**
- **NS5B, GPX, GPY** — corrected at **L2 (the router)**: it turns an unwarranted docking claim
  into an honest annotate_only / refuse / defer.
- **MPRO, PB1** — corrected at **L3 (the verification gate)**: it downgrades a docking ranking
  that does not beat the 2D-similarity baseline (or lacks a reproduced pose).

Run-to-run consistency: **1.0** (deterministic).

## Interpretation
- The raw pipeline (L0) emits a structure-based enrichment claim on **every** target — and on
  this held-out set that claim is warranted on **none** (each target either shouldn't be docked
  or loses to a trivial 2D baseline). That is precisely the epistemic failure Phase R prevents.
- **Dossier alone (L1) changes no decision** — an honest null: legibility without action does
  nothing. The value is in the router and the gate acting on the dossier.
- The router (L2) and the verification gate (L3) fix **disjoint** sets of targets; neither alone
  reaches full accuracy. This routing-type-dependence is the evidence that the integrating
  architecture — not the docking scorer — is the contribution.

## Honest limitations
- The ground-truth dispositions are **curated** (5 targets), not a large prospective set — this
  demonstrates the mechanism, it is not a population estimate. Ground truth for "was docking
  defensible here" is inherently sparse (the plan's stated caveat).
- The ablation measures **decision correctness**, deliberately not enrichment; it does not claim
  docking works, only that the workflow makes honest decisions about when it doesn't.
- R3 PlaybookMemory is still stubbed, so the "playbook prior" lever is not exercised here.
