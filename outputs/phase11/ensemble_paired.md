# Phase 11 WI-3b — Ensemble vs single-structure (paired bootstrap)

No significant ensemble effect on bedroc: paired Δ=0.0881, 95% CI [-0.172, 0.3708] includes 0 (P(Δ>0)=0.7666). Single-structure signal is stable under conformational sampling — established by a PAIRED test, not CI overlap.

| Metric | median Δ(ens−single) | 95% CI | P(Δ>0) | Holm reject |
|---|---|---|---|---|
| bedroc | 0.0881 | [-0.172, 0.3708] | 0.7666 | False |
| logauc | -0.0035 | [-0.0944, 0.0852] | 0.4661 | False |
| roc_auc | -0.1112 | [-0.2118, -0.0136] | 0.0135 | False |
| ef1 | 0.0 | [0.0, 2.1739] | 0.0891 | False |

This replaces the invalid CI-overlap comparison (Schenker & Gentleman 2001; Cumming 2009).
