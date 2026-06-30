# Research Accuracy Execution Plan

This plan turns the remaining accuracy-validation work into concrete execution stages.
The code seams are integrated; these steps verify real-tool behavior, quantify accuracy
impact, and prepare publication-grade runs.

## Goals

1. Confirm the integrated code matches the documented architecture.
2. Quantify each accuracy strategy's effect on control recovery and enrichment.
3. Update scoring and reporting defaults only when validation data supports it.
4. Prepare for publication-grade benchmarking, MD, and FEP/ABFE runs.

## Stage 1: Full Test Suite

Goal: confirm the current code passes in the local environment.

Workflow:

```bash
python -m pytest tests -q
```

Success criteria:

- All tests pass.
- No network or GPU dependency is required.
- Any failure is fixed before real-tool validation starts.

Status: completed locally on this workspace.

## Stage 2: Positive-Control Validation With Real Tools

Goal: verify known RdRp inhibitor recovery with real Vina and FPocket.

Workflow:

1. Install AutoDock Vina and FPocket.
2. Ensure `vina` and `fpocket` are discoverable on `PATH` or through local tool paths.
3. Run:

```bash
python scripts/validate_controls.py
```

Success criteria:

- Console output reports `VERDICT: PASS`.
- At least 3 of 4 known RdRp inhibitor controls recover in the expected top set.
- `outputs/validation_controls.json` reflects the current scoring path.

Status: pending real tool/environment confirmation.

## Stage 3: ADMET-AI Predictions

Goal: populate the HTML report with real ADMET predictions for top leads.

Workflow:

1. Install `admet_ai` in the active environment.
2. Run VTA-Agent on a reference genome.
3. Open `outputs/{run_id}_report.html`.
4. Confirm ADMET columns are populated.

Success criteria:

- Top leads include hERG, oral bioavailability, and solubility values.
- ADMET remains annotation-only unless a later validation pass justifies scoring changes.
- Report notes correctly describe ADMET interpretation.

Status: pending package installation and real run.

## Stage 4: DL Rescoring Impact

Goal: determine whether GNINA or Boltzina improves validation metrics.

Workflow:

1. Use a GPU host with GNINA and/or Boltzina installed.
2. Run the VTA pipeline on the reference target.
3. Preserve `cnn_score`, `cnn_affinity`, and `boltzina_score` annotations.
4. Test candidate score blends outside `rank.py`.
5. Re-run positive-control validation and enrichment metrics.
6. Commit ranking changes only if validation improves.

Success criteria:

- Baseline LE-led score remains the fallback.
- Any DL-weighted score improves control recovery or enrichment.
- Score changes are documented with before/after metrics.

Status: pending GPU host and tool installation.

## Stage 5: MD Smoke Test

Goal: confirm the MD validation path runs with short step counts.

Workflow:

1. Install OpenMM, OpenMMForceFields, OpenFF Toolkit, ParmEd, and MDAnalysis.
2. Set short-run controls, for example:

```bash
MD_STEPS=5000000 MD_EQUIL_STEPS=100000 vta run path/to/genome.fasta --include-md
```

3. Check `outputs/{run_id}/md/`.
4. Confirm the report renders MD verdicts or labelled skips.

Success criteria:

- MD files are created for selected candidates when dependencies are available.
- Failed candidates are labelled without crashing the graph.
- `md_validated_leads` appears in the final state/report.

Status: pending MD environment.

## Stage 6: Production MD

Goal: generate full MD validation data for top leads.

Workflow:

1. Select top leads from the reference run.
2. Run `vta run path/to/genome.fasta --include-md` with production step counts.
3. Inspect RMSD, persistent contacts, and MM-GBSA outputs.
4. Confirm the report's MD table is publication-ready.

Success criteria:

- Top leads have stable/unstable/moderate verdicts.
- MM-GBSA appears where topology and engine support it.
- MD reranking is documented and reproducible.

Status: pending GPU/long compute.

## Stage 7: FEP/ABFE Runner

Goal: prepare production absolute binding free energy calculations.

Workflow:

1. Install or provide an `openfe-abfe` compatible runner on `PATH`.
2. Verify the runner contract:

```text
openfe-abfe <protein_pdb> <ligand_smiles> <output_dir>
```

3. Ensure it emits JSON on stdout:

```json
{"delta_g": -8.4, "error": 0.7, "method": "OpenFE"}
```

4. Run:

```bash
vta run path/to/genome.fasta --include-fep
```

Success criteria:

- `fep_validated_leads` is populated.
- The report renders the FEP validation table.
- Skipped or failed FEP jobs are explicit.
- Compute cost and timeline are estimated before full lead-series execution.

Status: pending ABFE runner and compute.

## Update Rule

After each stage:

- Update `docs/research_accuracy_integration_progress.md`.
- Record exact command, environment, and result.
- Preserve output artifacts needed for reproducibility.
- Do not change scoring defaults without validation metrics.
