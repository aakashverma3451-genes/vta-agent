# Phase 10: Ensemble Docking — Mpro

Ensemble: 7L11:A (primary holo) + 6Y2E_apo, 7K3T_holo.
N compounds: 100. Excluded (all conformers failed): 0.
Box centers: per-conformer His41/Cys145 CA midpoint (crystal frames differ).

## Enrichment comparison

| Metric | Phase 9 (single) | Phase 10 (ensemble) | Delta | CI overlap? |
|--------|-----------------|---------------------|-------|-------------|
| BEDROC(α=20) | 0.6767 [0.357, 0.893] | 0.7848 [0.461, 0.932] | +0.1081 | True |
| ROC-AUC | 0.5796 [0.467, 0.690] | 0.4686 [0.347, 0.578] | -0.1110 | True |
| logAUC | 0.1994 | 0.2028 | +0.0034 | — |
| EF1% | 2.0 | 2.0 | +0.0000 | — |

## Ranking delta — top 10 Phase 9 actives

| Compound | Phase 9 rank | Ensemble rank | Delta | Improved? |
|----------|-------------|---------------|-------|-----------|
| EDG-MED-ba1ac7b9-11 | 1 | 6 | +5 | False |
| ALP-POS-ecbed2ba-1 | 2 | 4 | +2 | False |
| ALP-POS-6d96567b-2 | 4 | 10 | +6 | False |
| EDJ-MED-fa7708b3-3 | 5 | 34 | +29 | False |
| LUO-POS-868e8996-7 | 8 | 1 | -7 | True |
| ALP-POS-f1807566-1 | 9 | 5 | -4 | True |
| MAT-POS-e9e99895-3 | 11 | 37 | +26 | False |
| EDJ-MED-8bb691af-4 | 12 | 52 | +40 | False |
| MAT-POS-4223bc15-8 | 14 | 48 | +34 | False |
| MAT-POS-86c60949-2 | 17 | 11 | -6 | True |

## Interpretation

Ensemble (3 conformers) BEDROC 0.7848 vs Phase 9 single-structure 0.6767 (delta +0.1081). CI bands overlap — ensemble change is within noise; Phase 9 single-structure signal is stable w.r.t. conformational sampling. Top-10 Phase 9 actives: 3/10 improved rank in ensemble. Gate unchanged: DO NOT PROMOTE without leakage-controlled DL improvement at non-overlapping CIs.

Gate: **DO NOT PROMOTE** — no DL rescore wired; modest Mpro docking signal unchanged.
