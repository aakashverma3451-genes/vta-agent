# Phase 9C Gate Decision: SARS-CoV-2 Mpro

Decision: **DO NOT PROMOTE any DL or consensus ranking term.**

Evidence:

- The powered Mpro benchmark has not produced a number yet: real Vina receptor
  prep is blocked by a Meeko 0.7.1 Mpro-specific bug (see
  `outputs/phase9/mpro_noncovalent_benchmark.json`).
- With no powered enrichment result, there is no leakage-controlled DL-vs-baseline
  comparison to justify a ranking change.

Therefore consensus/DL outputs remain annotation-only and ranking weights are
unchanged. This negative result is logged deliberately: even a powered benchmark does
not, on its own, license a ranking change without a leakage-controlled DL improvement.
