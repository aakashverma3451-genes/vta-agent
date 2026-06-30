# Phase 8 Gate Decision: Matched Benchmark

Decision: DO NOT PROMOTE any DL or consensus term.

Evidence:

- Target: TiLV PB1 8PSO chain B.
- Benchmark file: `outputs/phase8/matched_benchmark_tilv_pb1.json`
- N = 54: 4 actives / 50 decoys.
- BEDROC(alpha=20) = 0.4062
- EF1% = 13.5
- logAUC = 0.3883
- ROC-AUC = 0.785

Reason:

- The benchmark is underpowered for ranking-term promotion.
- Decoys are presumed inactive, not experimentally verified inactive.
- GNINA/Boltz-2 training-target exclusion data were not supplied.
- No leakage-controlled DL-vs-baseline enrichment improvement has been shown.

Ranking remains unchanged.
