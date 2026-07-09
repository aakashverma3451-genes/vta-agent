# Phase S — activity-cliff benchmark (Moonshot Mpro)

Fair-arena test: on activity cliffs 2D-similarity is ~chance by design, so a structure-based win is possible and meaningful (van Tilborg 2022). NOT a Moonshot-whole enrichment claim (that is unwinnable by design).

- docked compounds with measured potency: **1109**
- activity-cliff pairs (Tanimoto ≥ 0.7, |ΔpIC50| ≥ 1.0): **1193** (powered)

## Cliff-pair ranking accuracy (median [95% CI]); chance = 0.50

| Method | accuracy | vs chance |
|---|---|---|
| **Vina (−ΔG)** | 0.4183 [0.3906, 0.4468] | significant |
| 2D-kNN pIC50 | 0.6655 [0.6387, 0.6924] | (near-chance by design) |

Paired Vina − 2D: Δ -0.2473, 95% CI [-0.2875, -0.207], P(Δ>0) 0.0.


## Strict-cliff variant (Tanimoto ≥ 0.9, 63 pairs)

Vina 0.4444 [0.3175, 0.5714] | 2D-kNN 0.8413 [0.746, 0.9206] | paired Vina−2D -0.3968 [-0.5556, -0.2222]

_Purer cliffs → 2D closer to chance; the Vina−2D gap is the fair-arena test at its strictest._

**Structure-based skill NOT demonstrated on activity cliffs: Vina does not beat both chance and the 2D-kNN baseline with non-crossing CIs on this powered cliff set.**

## Caveats
- Vina accuracy is significantly BELOW chance (0.42, CI upper bound < 0.5) — it is ANTI-correlated with potency on cliffs, ranking the less-potent analog as the stronger binder ~58% of the time. Most likely Vina's known size bias (larger, more elaborated but not more potent analogs get more-negative ΔG). Hypothesis, not proven.
- The 2D baseline is a kNN-pIC50 NEIGHBOURHOOD QSAR (excluding both pair members), NOT a pairwise-similarity test. It is therefore NOT forced to chance on cliffs — the strict Tanimoto≥0.9 variant made it BETTER (0.84), not worse, because Moonshot's dense congeneric series carry strong neighbourhood potency signal. So the honest claim is 'docking loses to a trivial ligand-based QSAR even on cliffs', not 'docking wins where similarity is disabled'. A pairwise-2D baseline would be ~0.5 by design.
- Cliff pairs are limited to compounds with BOTH a docked Vina ΔG AND a numeric measured IC50 (1,109 of the measured set); one consistent 7L11 real-Vina protocol.
- Measures cliff-pair RANKING, not enrichment — a null here does not claim docking is useless, only that it shows no marginal ranking skill in the fair arena.
