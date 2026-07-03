# Phase 11 WI-5 Job B — HCV NS5B NI, property-UNMATCHED decoys

Structure: 2XI3 chain A (NS5B + catalytic Mg2+). Ratio 638:43. Decoy recipe: ChEMBL drug-like (MW 250-500, ≤1 Ro5 violation), phosphate-free, scaffold-distinct — deliberately property-UNMATCHED to the charged triphosphate actives (Stein 2021 DUDE-Z philosophy). DeepCoy absent.

## Baselines vs Vina (median [95% CI])

| Method | BEDROC | logAUC | ROC-AUC | EF1% |
|---|---|---|---|---|
| random | 0.0847 [0.0244, 0.1722] | 0.144 [0.1044, 0.1931] | 0.5006 [0.4114, 0.5881] | 0.0 [0.0, 4.5249] |
| 2d_sim | 1.0 [1.0, 1.0] | 1.0 [1.0, 1.0] | 1.0 [1.0, 1.0] | 15.8372 [12.3818, 21.9677] |
| vina | 0.0807 [0.0323, 0.1532] | 0.2305 [0.2004, 0.2651] | 0.7675 [0.7216, 0.8095] | 0.0 [0.0, 0.0] |

Paired Vina−2Dsim (BEDROC): Δ -0.9197, 95% CI [-0.9671, -0.8475].

## Interpretation

This benchmark DOES run with a property-unmatched recipe (unlike the property-matched recipe, which could not build a decoy set). Read the numbers with the DUDE-Z caveat: high enrichment against neutral drug-like decoys may reflect charge/size discrimination of the triphosphate actives rather than pocket-specific binding — the paired Vina-vs-2D-similarity test separates those. Reported as-is.
