# Phase S — activity-cliff benchmark: DONE (result: directional sign-flip, underpowered)

**Date:** 2026-07-07 · **Owner:** session-a (opus)

**Goal:** test docking on the one fair arena where a structure-based win is possible — activity
cliffs (near-identical in 2D, large ΔpIC50), where 2D-similarity is ~chance by design. NOT an
attempt to beat 2D on the analog-clustered Moonshot set as a whole (unwinnable; Boby 2023).

## Result (`outputs/phaseS/activity_cliffs.{json,md}`)
Built from committed data only: Job A 31:1 real-Vina ΔG (764 docked compounds) joined to the
stratified dataset's measured IC50 + SMILES. Cliff pairs = ECFP4 Tanimoto ≥ 0.7 AND |ΔpIC50| ≥ 1.

| Method | cliff-ranking accuracy (median [95% CI]) | vs chance (0.50) |
|---|---|---|
| **Vina (−ΔG)** | **0.647 [0.412, 0.882]** | above chance, but CI crosses 0.5 |
| 2D-kNN pIC50 | 0.353 [0.118, 0.588] | **below** chance (cliff pathology) |

Paired **Vina − 2D = +0.294**, 95% CI **[−0.118, 0.706]**, P(Δ>0) = 0.92. **n = 17 cliff pairs
(UNDERPOWERED).**

## Verdict — honest
- **Not demonstrated:** every CI crosses its null, so no claim clears significance. Reported as
  the honest null (`structure_based_skill_on_cliffs_demonstrated = false`).
- **But the sign flips in the fair arena.** On enrichment, 2D crushed Vina (Job A paired Δ
  BEDROC −0.24, P=0.01). On cliffs, **Vina points *above* chance (0.65) and 2D points *below*
  (0.35)**, with Vina out-ranking 2D by ~29 points (P(Δ>0)=0.92). This is the first time in the
  project that docking has pointed the *right* way relative to 2D — consistent with the
  hypothesis that docking's value, if any, is in congeneric/cliff ranking, not enrichment
  (van Tilborg 2022). It is a **hypothesis-generating signal, not a result.**

## Why underpowered (expected, pre-registered as ~likely)
Moonshot is analog-clustered but the *docked* subset with a numeric measured IC50 yields only
17 pairs meeting the cliff criteria — the power is bounded by the docked pool, exactly the
"demonstration-grade / CIs wide" outcome flagged before running (cf. TiLV).

## What this changes
- **Evidentiary (guaranteed):** the last "unfair benchmark" objection is closed — there is now
  a fair-arena test in the ladder, and a reusable cliff-scoring module (`vta/eval/cliffs.py`)
  that feeds the same paired honesty gate.
- **A discriminating yardstick now exists** for a learned rescorer: run GNINA/RTMScore on these
  exact 17 pairs (+ an augmented set) and see whether the directional signal becomes significant.
- **No ranking change** — the gate keeps docking annotation-grade; nothing here promotes it.

## Next steps (not taken this session)
1. **Power the cliff set:** augment with ChEMBL Mpro congeneric series (more docked pairs with
   measured IC50) so the ±0.29 Vina−2D gap can be tested at significance.
2. **Rescorer on the cliff set:** GNINA CNNaffinity / RTMScore vs the same 2D-kNN baseline on the
   identical pairs — the honest path to a *powered* structure-based positive.
