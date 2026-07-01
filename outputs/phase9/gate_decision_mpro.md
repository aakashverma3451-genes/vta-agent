# Phase 9/10 Gate Decision: SARS-CoV-2 Mpro

Decision: **DO NOT PROMOTE any DL or consensus ranking term.**

## Phase 9 single-structure (7L11:A)

- BEDROC 0.68 [0.357, 0.893], ROC-AUC 0.58 [0.467, 0.690]
- ROC-AUC CI crosses 0.5 → modest signal, not a promotion basis.

## Phase 10 ensemble (7L11:A + 6Y2E:A apo + 7K3T:A holo)

- BEDROC 0.7848 [0.461, 0.932], ROC-AUC 0.4686 [0.347, 0.578]
- BEDROC delta +0.1081, ROC-AUC delta -0.1110
- 3/10 top Phase 9 actives improved rank in ensemble
- Phase 10 ensemble BEDROC CI [0.461, 0.932] overlaps Phase 9 BEDROC CI [0.357, 0.893] — conformational sampling does not significantly change enrichment. Phase 9 single-structure signal is stable w.r.t. the tested ensemble.

## Promotion criteria (unchanged)

- A powered enrichment result exists (Phase 9), but promotion requires a DL/consensus
  signal that beats the Vina baseline with NON-OVERLAPPING CIs.
- `consensus_node` reads `cnn_affinity`/`boltzina_score`, which are not populated
  (no DL rescore is wired), so no DL term can be evaluated against the baseline.
- Ensemble docking (Phase 10) quantifies conformational stability — it does not
  itself constitute a DL improvement.

Therefore consensus/DL outputs remain annotation-only and ranking weights are
unchanged. This negative result is logged deliberately.
