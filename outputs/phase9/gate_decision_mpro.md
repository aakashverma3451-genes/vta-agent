# Phase 9C Gate Decision: SARS-CoV-2 Mpro

Decision: **DO NOT PROMOTE any DL or consensus ranking term.**

Evidence:

- A powered Mpro enrichment result exists, but promotion requires a DL/consensus
  signal that beats the Vina baseline with NON-OVERLAPPING CIs.
- `consensus_node` reads `cnn_affinity`/`boltzina_score`, which are not populated
  (no DL rescore is wired), so no DL term can be evaluated against the baseline.

Therefore consensus/DL outputs remain annotation-only and ranking weights are
unchanged. This negative result is logged deliberately: even a powered benchmark does
not, on its own, license a ranking change without a leakage-controlled DL improvement.
