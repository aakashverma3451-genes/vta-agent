# Phase 8 Activation Report

## Target Decision

Phase 8 is activated on one validation target:

- Virus/target: TiLV PB1
- Structure: 8PSO chain B
- Pocket: NTP catalytic site
- Rationale: this is the existing validated real-tool path for FPocket and AutoDock Vina.

## 8A Benchmark Evidence

Two benchmark artifacts now exist:

1. `outputs/phase8/benchmark_tilv_pb1.json`
   - Uses the pre-existing real Vina/FPocket control run and mechanism-distinct demo decoys.
   - Verdict: demonstration only; property-matched decoys still required.

2. `outputs/phase8/matched_benchmark_tilv_pb1.json`
   - Uses real ChEMBL API property-matched presumed decoys fetched by
     `scripts/fetch_phase8_decoys.py`.
   - Scores controls plus matched decoys through the real TiLV PB1 Vina/FPocket path using
     `scripts/run_phase8_matched_benchmark.py`.
   - N = 54: 4 actives and 50 decoys.
   - BEDROC(alpha=20) = 0.4062
   - EF1% = 13.5
   - logAUC = 0.3883
   - ROC-AUC = 0.785
   - Verdict: defensible demonstration; still underpowered for publication.

The matched decoy set is real but not publication-grade: it contains 50 presumed
decoys total, not the requested 30-50 decoys per active, and the decoys are not
experimentally verified inactives.

## 8B Locked Predictions

Locked prediction registry:

- File: `outputs/phase8/locked_predictions_tilv_pb1.json`
- Run ID: `phase8_tilv_pb1`
- Locked at: 2026-06-29T09:19:32.396246+00:00
- Top locked prediction: Ribavirin, score 0.635, input hash `38f398344ae130ee`

Experimental plan:

- File: `outputs/phase8/experimental_plan_tilv_pb1.md`
- Primary assay: biochemical TiLV PB1/RdRp inhibition IC50 if available.
- Alternative assay: cell EC50 plus CC50 to compute selectivity index.
- Success criterion: at least one locked top prediction shows reproducible dose-response;
  for cell assays, CC50/EC50 should exceed 10.

## 8E Gate Decision

No DL term is promoted.

Reasons:

- The matched benchmark is still underpowered.
- GNINA/Boltz-2 scorer training-target exclusion data were not supplied.
- Consensus/DL annotations did not undergo a leakage-controlled enrichment comparison
  against this matched benchmark.

Ranking remains unchanged.

## Honest Control-Gate Before/After

- Original real-tool control gate: 3/4 controls in top 5, PASS.
- Phase 8 expanded matched-decoy benchmark: 1/4 controls in top 5.

The drop is recorded as evidence, not tuned away.

## Still Skipped

- Full property-matched benchmark: needs 30-50 docked decoys per active.
- Selectivity panel: needs Foldseek/HHpred human homolog panel and anti-target structures.
- Resistance panel: needs curated mutation positions or entropy/dN/dS-derived mutants.
- ADMETlab 3.0: provider activation still needs API integration and live calls.
- KG/literature slice: needs PrimeKG/DRKG source data and Europe PMC/PubTator runs.
