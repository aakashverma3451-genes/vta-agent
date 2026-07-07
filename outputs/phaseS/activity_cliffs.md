# Phase S — activity-cliff benchmark (Moonshot Mpro)

Fair-arena test: on activity cliffs 2D-similarity is ~chance by design, so a structure-based win is possible and meaningful (van Tilborg 2022). NOT a Moonshot-whole enrichment claim (that is unwinnable by design).

- docked compounds with measured potency: **764**
- activity-cliff pairs (Tanimoto ≥ 0.7, |ΔpIC50| ≥ 1.0): **17** (UNDERPOWERED — demonstration-grade)

## Cliff-pair ranking accuracy (median [95% CI]); chance = 0.50

| Method | accuracy | vs chance |
|---|---|---|
| **Vina (−ΔG)** | 0.6471 [0.4118, 0.8824] | CI crosses 0.5 |
| 2D-kNN pIC50 | 0.3529 [0.1176, 0.5882] | (near-chance by design) |

Paired Vina − 2D: Δ 0.2941, 95% CI [-0.1176, 0.7059], P(Δ>0) 0.8866.

**Structure-based skill NOT demonstrated on activity cliffs: Vina does not beat both chance and the 2D-kNN baseline with non-crossing CIs — and with only 17 cliff pairs the test is underpowered (demonstration-grade, CIs wide).**

## Caveats
- Cliff pairs are limited to compounds that were BOTH docked (have Vina ΔG) AND carry a numeric measured IC50 — a subset of Moonshot, so power is bounded by the docked set.
- The 2D baseline is a kNN-pIC50 predictor excluding both pair members (its fairest estimate); by construction it is ~chance on true cliffs — that is the point.
- Measures cliff-pair RANKING, not enrichment; a null here does not claim docking is useless, only that it shows no marginal ranking skill even where 2D is disabled.
