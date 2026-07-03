# Phase 11 WI-4 Job A — Mpro non-covalent at 740:24

Real Vina: True. 24 actives / 740 measured inactives. EF1% ceiling now 31.8.

## Baselines vs Vina (median [95% CI])

| Method | BEDROC | logAUC | ROC-AUC | EF1% |
|---|---|---|---|---|
| random | 0.0625 [0.007, 0.1598] | 0.1436 [0.0928, 0.2099] | 0.4997 [0.3842, 0.6189] | 0.0 [0.0, 7.9583] |
| 2d_sim | 0.398 [0.2369, 0.5568] | 0.3835 [0.2891, 0.4948] | 0.8197 [0.7221, 0.8997] | 9.55 [0.0, 21.7045] |
| vina | 0.153 [0.0568, 0.2779] | 0.2213 [0.151, 0.3029] | 0.634 [0.5066, 0.754] | 0.0 [0.0, 9.55] |

Paired Vina−2Dsim (BEDROC): Δ -0.2446, 95% CI [-0.4385, -0.0444], P(Δ>0) 0.0099.

**At 30:1, Vina still does NOT beat 2D-similarity on BEDROC — structure-based enrichment not demonstrated even at a non-degenerate ratio.**
