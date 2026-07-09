# Phase S — activity-cliff benchmark: DONE & POWERED (result: docking loses to trivial QSAR)

**Date:** 2026-07-07 · **Owner:** session-a (opus)

**Goal:** test docking on activity cliffs — the fair arena where 2D-similarity reasoning is
weakest — on a POWERED set, after docking 345 more compounds to lift the cliff pool from 17 to
1,193 pairs.

## Powered result (`outputs/phaseS/activity_cliffs.{json,md}`)
1,109 docked compounds with measured IC50 → **1,193 cliff pairs** (Tanimoto ≥ 0.7, |ΔpIC50| ≥ 1).
All 345 new docks succeeded (0 fail, 0 timeout), one consistent 7L11 real-Vina protocol.

| Method | cliff-ranking accuracy (median [95% CI]); chance = 0.50 |
|---|---|
| **Vina (−ΔG)** | **0.418 [0.391, 0.447]** — significantly **BELOW** chance |
| 2D-kNN pIC50 | **0.666 [0.639, 0.692]** — significantly above chance |

Paired **Vina − 2D = −0.247**, 95% CI **[−0.288, −0.207]**, P(Δ>0) = 0. **Powered (n=1,193).**

**Strict variant (Tanimoto ≥ 0.9, 63 pairs):** Vina 0.44 [0.32, 0.57], 2D-kNN **0.84** [0.75, 0.92],
paired −0.40 [−0.56, −0.22].

## Verdict — powered, and it reverses the earlier fluke
- **Structure-based skill NOT demonstrated — and worse than that:** Vina is *significantly below
  chance* (0.42) on cliff-pair ranking, i.e. **anti-correlated with potency** — it ranks the
  less-potent analog as the stronger binder ~58% of the time. The trivial 2D-kNN QSAR beats it
  by ~25 points with a tight CI.
- **This overturns the 17-pair pilot** (which showed Vina 0.65 > 2D 0.35, paired +0.29). That
  was a small-sample fluke; powering the benchmark reversed the sign decisively. **This is
  exactly the value of powering** — it caught and corrected a misleading positive.

## Two honest scientific points (both self-corrections)
1. **Why Vina is *below* chance:** most likely **Vina's known size bias** — larger, more
   elaborated analogs get more-negative ΔG regardless of true potency, so within a congeneric
   pair Vina systematically favours the wrong member. Hypothesis, not proven here.
2. **The pre-registered expectation was WRONG:** I expected stricter cliffs (Tanimoto ≥ 0.9) to
   push the 2D baseline toward chance. The opposite happened — 2D-kNN got **better (0.84)**. The
   reason: the kNN baseline is a **neighbourhood QSAR** (excluding both pair members), not a
   pairwise-similarity test, so on Moonshot's dense congeneric series it exploits strong
   neighbourhood potency signal and is NOT disabled by cliffs. So the honest headline is
   **"docking loses to a trivial ligand-based QSAR even on activity cliffs,"** not "docking wins
   where similarity is disabled." A pairwise-2D baseline would be ~0.5 by construction.

## What this changes
- The project's core negative is now **powered in the fair arena**: docking shows no
  structure-based ranking skill on cliffs — it is beaten by trivial QSAR and is itself
  anti-correlated with potency. Stronger and more defensible than the enrichment result alone.
- **No ranking change** — the R4 gate keeps docking annotation-grade; nothing here promotes it.
- The cliff benchmark + 1,109-compound docked cache are now a reusable, powered **yardstick for a
  rescorer**: GNINA/RTMScore must beat *both* this 2D-kNN QSAR and chance on these exact pairs.

## Next step (the honest path to any structure-based positive)
Run **GNINA CNNaffinity / RTMScore** on these 1,193 cliff pairs vs the same 2D-kNN baseline
through the paired gate. If a learned rescorer corrects Vina's size-bias anti-correlation and
beats the QSAR, that is a real, powered, first structure-based win. If not, the honest ceiling
stands.
