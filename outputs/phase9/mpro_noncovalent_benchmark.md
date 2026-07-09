# Phase 9C Benchmark: SARS-CoV-2 Mpro (non-covalent)

Target: SARS-CoV-2 Mpro (3CLpro), non-covalent. Structure: 7L11 chain A.
Binding mode: non_covalent actives only (covalent held out).
Control set: experimentally measured Moonshot inactives (real, not presumed decoys).
Engine: AutoDock Vina (real Vina: True; receptor prep: openbabel_fallback).
Verdict: powered benchmark — real measured inactives, informative CIs (not [0,1]); docking signal MODEST — ROC-AUC CI [0.4669, 0.6903] includes 0.5, EF1% ceiling-limited at 50/50.

## Primary metrics (point + bootstrap CI)

- N: 100 (50 actives / 50 inactives)
- BEDROC(alpha=20): 0.6767 — median 0.6823, 95% CI [0.3574, 0.8933]
- EF1%: 2.0 — median 1.9608, 95% CI [0.0, 2.439]
- logAUC: 0.1994 — median 0.2056, 95% CI [0.137, 0.3012]
- ROC-AUC: 0.5796 — median 0.5795, 95% CI [0.4669, 0.6903] (secondary)
- CIs informative (not [0,1]): True

## Control-draw spread (robustness to inactive subsample)

- bedroc: median 0.7177 (min 0.6889, max 0.8154)
- EF1%: median 1.8 (min 1.8, max 1.8)
- log_auc: median 0.1976 (min 0.1916, max 0.2283)
- roc_auc: median 0.5795 (min 0.5745, max 0.593)

## Leakage audit

- Status: pass_no_overlap
- Overlap: []

Ranking weights unchanged. Consensus/DL outputs remain annotation-only.
