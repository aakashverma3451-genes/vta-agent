# Phase 9C Benchmark: SARS-CoV-2 Mpro (non-covalent) — BLOCKED

Status: **blocked_receptor_prep**.
Verdict: blocked — powered benchmark not yet runnable; no number claimed.

## Blocker

AutoDock Vina is installed, but Meeko 0.7.1 receptor preparation fails deterministically on the Mpro chain (mk_prepare_receptor raises 'update_H_positions: Updated 1 H positions but deleted N' for every Mpro PDB tried: 6Y2E, 6M03, 7K3T, 7TLL, 7L10, 6W63, 7RFW, 7VH8, 7L11, 5R8T, 7BB2). The same code path preps the TiLV 8PSO receptor successfully, so the engine and pipeline are intact; the failure is Meeko-vs-Mpro specific. No fallback receptor prep tool (OpenBabel / reduce / ADFR) is installed.

## What is ready (no docking required)

- Dataset: `vta/data/mpro/mpro_dataset.json`
- Stratified (binding mode): `vta/data/mpro/mpro_stratified.json`
- Selected non-covalent actives: 50
- Selected measured inactives: 50
- Active-site wiring: EXPERIMENTAL_PDB/EXPERIMENTAL_ACTIVE_SITE += MPRO (7L11:A)

## Next step

Provide a working receptor-prep path (OpenBabel fallback seam, or a Meeko version without this bug), then re-run scripts/run_mpro_benchmark.py to produce the powered, CI-bearing enrichment number. No mock/fabricated number is emitted.

No mock or fabricated enrichment number is reported for Mpro.
