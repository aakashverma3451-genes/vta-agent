# Phase 11 — Trivial baselines, paired significance, LE gate

Every number is a bootstrap median + 95% CI (n=10000). Significance is by paired bootstrap of the difference on the SAME split — never CI overlap.

## Mpro_noncovalent (powered; 50 actives / 50 inactives)

### WI-1 baselines vs Vina (median [95% CI])

| Method | BEDROC | logAUC | ROC-AUC | EF1% |
|---|---|---|---|---|
| random | 0.5007 [0.2059, 0.7916] | 0.1494 [0.0978, 0.2315] | 0.4996 [0.3864, 0.6128] | 0.0 [0.0, 2.0] |
| 2d_sim | 0.9171 [0.7229, 0.9868] | 0.3621 [0.2596, 0.4994] | 0.7572 [0.6554, 0.8463] | 2.0 [1.6667, 2.439] |
| vina | 0.6831 [0.3514, 0.9054] | 0.2057 [0.1348, 0.3077] | 0.5808 [0.4667, 0.6923] | 1.9608 [0.0, 2.439] |

### WI-3 paired test — Vina vs 2D-similarity (primary: BEDROC)

- median Δ(Vina−2Dsim) = -0.2238, 95% CI [-0.568, 0.0555], P(Δ>0) = 0.0594
- **structure-based enrichment NOT demonstrated on this target: Vina does not beat 2D-similarity memorization on BEDROC with a positive paired-difference CI.**

### WI-6 LE gate (primary: BEDROC)

- composite vs ΔG-only: median Δ = -0.2653, CI [-0.5791, 0.0581], P(Δ>0) = 0.0523
- **Gate: DEMOTE LE to secondary annotation; make ΔG the primary ranking term** (ranking change applied: True)

## TiLV_PB1 (UNDERPOWERED; 4 actives / 50 inactives)

### WI-1 baselines vs Vina (median [95% CI])

| Method | BEDROC | logAUC | ROC-AUC | EF1% |
|---|---|---|---|---|
| random | 0.0277 [0.0, 0.4441] | 0.1351 [0.0378, 0.3784] | 0.5 [0.21, 0.79] | 0.0 [0.0, 13.5] |
| 2d_sim | 1.0 [1.0, 1.0] | 1.0 [1.0, 1.0] | 1.0 [1.0, 1.0] | 13.5 [6.75, 54.0] |
| vina | 0.0974 [0.0, 0.4793] | 0.1735 [0.0146, 0.4157] | 0.5109 [0.0962, 0.9281] | 0.0 [0.0, 0.0] |

### WI-3 paired test — Vina vs 2D-similarity (primary: BEDROC)

- median Δ(Vina−2Dsim) = -0.9026, 95% CI [-1.0, -0.5117], P(Δ>0) = 0.0
- **structure-based enrichment NOT demonstrated on this target: Vina does not beat 2D-similarity memorization on BEDROC with a positive paired-difference CI.**

### WI-6 LE gate (primary: BEDROC)

- composite vs ΔG-only: median Δ = 0.2689, CI [-0.3802, 1.0], P(Δ>0) = 0.6632
- **Gate: DEMOTE LE to secondary annotation; make ΔG the primary ranking term** (ranking change applied: True)

