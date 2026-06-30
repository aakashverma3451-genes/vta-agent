# Claude Code Handoff — VTA-Agent Scientific Validation

Date: 2026-06-29  
Workspace: `/Users/dna/Desktop/vtaagent/vta-agent`

## Current Position

VTA-Agent has moved from architecture work into evidence generation. Phases 0-8 are
implemented as guarded scientific seams, and Phase 9 has started.

The most important result so far is not a high benchmark number. It is that the pipeline
successfully deflated its own optimistic validation:

- Original 9-ligand real-tool control gate: 3/4 controls in top 5.
- Phase 8 matched-decoy TiLV PB1 benchmark: 1/4 controls in top 5.
- Bootstrap CI on TiLV PB1 BEDROC: `[0.0009, 1.0]`.

Interpretation: TiLV PB1 with 4 actives is an underpowered demonstration, not a claim-grade
benchmark. Do not optimize or tune to recover the old 3/4 number.

## Non-Negotiable Constraints

- Do not change `vta/nodes/rank.py` weights unless a leakage-controlled benchmark with
  non-overlapping CIs justifies it.
- Do not pool HCV NS5B NI and NNI actives into one pocket benchmark.
- Do not fabricate triphosphate active-form SMILES, metal coordinates, decoys, or assay data.
- Missing data/tool means labelled skip with provenance.
- Existing `.claude/COLLAB.md` deletion is unrelated; leave it alone.
- Keep tests green.

## Test Status

Last full suite:

```bash
PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python -m pytest tests -q
```

Result:

```text
169 passed in 27.24s
```

## Key Artifacts

### Phase 8 TiLV PB1 Evidence

- `outputs/phase8/matched_benchmark_tilv_pb1.json`
- `outputs/phase8/benchmark_matched_tilv_pb1.md`
- `outputs/phase8/benchmark_ci_tilv_pb1.md`
- `outputs/phase8/gate_decision_matched_tilv_pb1.md`
- `outputs/phase8/locked_predictions_tilv_pb1.json`
- `docs/phase8_activation_report.md`

Current TiLV PB1 matched benchmark:

- Target: TiLV PB1, 8PSO chain B.
- N = 54: 4 actives / 50 ChEMBL property-matched presumed decoys.
- BEDROC(alpha=20): 0.4062.
- EF1%: 13.5.
- logAUC: 0.3883.
- ROC-AUC: 0.785.
- Controls in top 5: 1/4.

Bootstrap CIs:

- BEDROC median 0.4062, 95% CI `[0.0009, 1.0]`.
- EF1% median 9.0, 95% CI `[0.0, 39.825]`.
- logAUC median 0.3877, 95% CI `[0.1327, 1.0]`.
- ROC-AUC median 0.7921, 95% CI `[0.5976, 1.0]`.

Conclusion: underpowered demonstration only.

### Phase 9 HCV NS5B Actives

- `vta/data/phase9_actives/hcv_ns5b.json`
- `vta/data/phase9_actives/hcv_ns5b_stratified.json`
- `vta/data/phase9_actives/hcv_ns5b_ni_active_site.json`
- `vta/data/phase9_actives/hcv_rdrp_mechanism_target.json`
- `vta/data/phase9_actives/influenza_pa_pb1.json`
- `docs/phase9_benchmark_maturity.md`

Important HCV stratification result:

Raw `CHEMBL5375` is not a clean same-pocket active set.

Counts from `hcv_ns5b_stratified.json`:

- `NI`: 1
- `NNI`: 6
- `NS5B_site_unknown`: 8
- `indirect_antiviral_assay`: 35
- `off_target_hcv_drug`: 2

Correct active-site scope:

- File: `vta/data/phase9_actives/hcv_ns5b_ni_active_site.json`
- Source: ChEMBL `CHEMBL4296320`
- Class: NI / nucleoside-like active-site inhibitors
- Count: 43 actives
- Pocket: NS5B catalytic active site

This is the next benchmark input. Do not use raw `CHEMBL5375` pooled against one pocket.

## Code Added Recently

Phase 8:

- `scripts/phase8_activate.py`
- `scripts/fetch_phase8_decoys.py`
- `scripts/run_phase8_matched_benchmark.py`
- `tests/test_phase8_activation.py`

Phase 9:

- `scripts/fetch_phase9_actives.py`
- `scripts/stratify_hcv_ns5b_actives.py`
- `scripts/build_hcv_ns5b_ni_scope.py`
- `tests/test_phase9_actives.py`
- Bootstrap CI support in `vta/eval/metrics.py`

Earlier architecture seams:

- `vta/eval/decoys.py`
- `vta/eval/leakage.py`
- `vta/eval/splits.py`
- `vta/nodes/consensus.py`
- `vta/nodes/selectivity.py`
- `vta/nodes/resistance.py`
- `vta/nodes/target_prioritization.py`
- `vta/nodes/structure_qc.py`
- `vta/experiment/registry.py`
- `vta/kg/schema.py`
- `vta/lit/extraction.py`

## Scientific Caveats To Preserve

### TiLV PB1

TiLV PB1 still docks parent forms for remdesivir, sofosbuvir, and molnupiravir because
curated triphosphate active-form SMILES are not committed. The 1/4 result may reflect both
scoring limitations and mechanistic-input limitations.

### HCV NS5B

NS5B inhibitors bind different pockets:

- NI/nucleotide/nucleoside inhibitors bind the catalytic active site.
- NNIs bind allosteric thumb/palm pockets.

Therefore:

- Headline HCV benchmark should use `hcv_ns5b_ni_active_site.json`.
- Dock those actives to the NS5B catalytic active site.
- Optional NNI benchmark must be separate and allosteric-pocket-specific.

### Metrics

Report intervals, not just points. With low actives, bootstrap CIs are themselves unstable,
but that instability is the result: the benchmark cannot estimate enrichment.

## Next Work — Phase 9B-3

Goal: build property-matched, scaffold-distinct decoys for the 43 HCV NS5B NI active-site
actives.

Input:

```text
vta/data/phase9_actives/hcv_ns5b_ni_active_site.json
```

Recommended output:

```text
vta/data/decoys_cache/hcv_ns5b_ni_matched.smi
outputs/phase9/hcv_ns5b_ni_decoy_provenance.json
outputs/phase9/hcv_ns5b_ni_matching_quality.json
```

Minimum target:

- Ideally 30 decoys per active: 43 * 30 = 1290 decoys.
- If ChEMBL API limits make that impractical, produce a smaller labelled demonstration set
  first, but do not call it publication-grade.

Matching properties:

- MW
- logP
- HBD
- HBA
- charge if available
- rotatable bonds
- scaffold distinct from active set

Useful existing code to adapt:

```text
scripts/fetch_phase8_decoys.py
vta/eval/splits.py
vta/eval/decoys.py
```

## Next Work — Phase 9B-4

After decoys exist, run the HCV NS5B NI active-site benchmark.

Required missing pieces:

1. Pick an HCV NS5B active-site structure from PDB.
2. Configure the active-site pocket.
3. Use active-form/triphosphate species where curated. Do not fabricate missing forms.
4. Dock NI actives and matched decoys through Vina.
5. Compute metrics with bootstrap CIs.

Recommended outputs:

```text
outputs/phase9/hcv_ns5b_ni_benchmark.json
outputs/phase9/hcv_ns5b_ni_benchmark.md
outputs/phase9/two_target_ci_table.md
outputs/phase9/gate_decision_hcv_ns5b_ni.md
```

The two-target table should compare:

- TiLV PB1: labelled underpowered demonstration.
- HCV NS5B NI active-site: powered benchmark, if enough decoys and valid structure/pocket.

## Commands

Run tests:

```bash
PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python -m pytest tests -q
```

Regenerate TiLV Phase 8 summary:

```bash
PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python scripts/phase8_activate.py
```

Fetch HCV NS5B actives:

```bash
PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python scripts/fetch_phase9_actives.py \
  --target-id CHEMBL5375 \
  --query 'Hepatitis C virus NS5B RNA-dependent RNA polymerase' \
  --label hcv_ns5b \
  --pages 5
```

Fetch HCV RNA polymerase mechanism-target actives:

```bash
PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python scripts/fetch_phase9_actives.py \
  --target-id CHEMBL4296320 \
  --query 'RNA-directed RNA polymerase inhibitor mechanism target' \
  --label hcv_rdrp_mechanism_target \
  --pages 5
```

Stratify raw HCV NS5B actives:

```bash
PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python scripts/stratify_hcv_ns5b_actives.py
```

Build HCV NS5B NI active-site scope:

```bash
PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python scripts/build_hcv_ns5b_ni_scope.py
```

## What Not To Do Next

- Do not move to Phase 10 before generating HCV NI decoys and at least a labelled benchmark
  attempt.
- Do not promote consensus/GNINA/Boltzina.
- Do not claim TiLV PB1 enrichment.
- Do not use all 52 CHEMBL5375 rows as active-site actives.
- Do not hide the 3/4 to 1/4 collapse; it is the strongest methodological evidence so far.

## Recommended Next Claude Code Prompt

```text
Continue VTA-Agent Phase 9B from docs/claude_code_handoff_2026-06-29.md.

Start with 9B-3:
- Build property-matched, scaffold-distinct decoys for
  vta/data/phase9_actives/hcv_ns5b_ni_active_site.json.
- Write decoys to vta/data/decoys_cache/hcv_ns5b_ni_matched.smi.
- Write provenance and matching-quality reports under outputs/phase9/.
- Do not fabricate decoys; use ChEMBL/DeepCoy/DEKOIS-style matching and label any
  underpowered result.
- Add hermetic tests for matching and provenance.
- Run PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python -m pytest tests -q.

Do not run the HCV benchmark until the decoy set and matching-quality report exist.
```
