# Research Accuracy Integration Progress

This document tracks progress against the accuracy-improvement ideas referenced by
`INTEGRATION_PLAN.md` and the missing source document named
`VTA-Agent - Research Papers That Can Improve Accuracy.md`.

That exact research-paper file is not present in this repository. This status is based
on the integration plan, current code, tests, and architecture documents.

Operational execution steps for real-tool validation are tracked in
`docs/research_accuracy_execution_plan.md`.

## Current Summary

Most of the accuracy architecture has been integrated as code-level seams. The default
pipeline can run without heavy tools, while GPU or external-tool features activate when
their dependencies exist.

Phase 0 safety rails are implemented: score provenance metadata, chemistry reality
flags, and active-species/prodrug annotations are added to lead records and reports.
Phase 1 scaffolding is also implemented: active-species resolution now runs upstream
of docking, docking records carry the selected docking species, and metal/ensemble
docking seams are labelled. Ranking remains unchanged.

Current default graph:

```text
classify -> router -> structure -> proteinttt -> pockets -> conservation
-> species_resolution -> dock -> conservation_contacts -> rescore -> boltzina
-> rank -> chemistry -> admet -> report
```

Optional MD graph:

```text
... -> admet -> md_select -> md_simulate -> md_analyze -> md_rerank -> report
```

Optional FEP graph:

```text
... -> md_rerank -> fep -> report
```

## Progress Table

| Accuracy Area | Implementation Status | Code | Notes |
|---------------|-----------------------|------|-------|
| ADMET-AI safety/drug-likeness | Integrated, annotation-only | `vta/nodes/admet.py` | Adds hERG, oral bioavailability, solubility when `admet_ai` is installed. |
| P2Rank pocket consensus | Integrated seam | `vta/nodes/pockets.py` | Adds `consensus` and `detectors` when P2Rank is available. |
| GNINA CNN rescoring | Integrated seam, annotation-only | `vta/nodes/rescore.py` | Requires `GNINA_BIN` or `gnina`; does not affect ranking yet. |
| Boltzina DL affinity | Integrated seam, annotation-only | `vta/nodes/boltzina.py` | Package API still needs verification against the real released package. |
| Boltz-2 structure fallback | Integrated seam | `vta/nodes/structure.py` | Used for large orphan proteins when `boltz` is available. |
| ProteinTTT fold refinement | Integrated seam | `vta/nodes/proteinttt.py` | Refines only low-pLDDT ESMFold structures. |
| Conservation JSD | Integrated | `vta/nodes/conservation.py` | Uses committed homolog MSA when available; otherwise labelled neutral fallback. |
| Ligand-contact conservation | Integrated | `vta/nodes/conservation_contacts.py` | Uses real pose contacts to make conservation ligand-specific. |
| MD validation | Integrated opt-in seam | `vta/nodes/md_*.py` | Requires OpenMM/MDAnalysis for real production runs. |
| MM-GBSA | Integrated analysis seam | `vta/nodes/md_analyze.py` | Runs when MM-GBSA engine and topology are available. |
| Enrichment metrics | Integrated | `vta/eval/metrics.py` | Provides EF, BEDROC, ROC-AUC. |
| Benchmark driver | Integrated | `scripts/benchmark_enrichment.py` | Uses current scores and demo decoys; real DUD-E/LIT-PCBA benchmark data still needed. |
| Positive-control validation | Integrated | `scripts/validate_controls.py` | Checks known RdRp inhibitor recovery. |
| FEP / ABFE | Integrated opt-in seam | `vta/nodes/fep.py` | Requires an `openfe-abfe` runner on PATH; production calculations still require expert setup and compute. |
| Chemistry filters / active species | Integrated, annotation-only | `vta/chem/*.py`, `vta/nodes/chemistry.py` | Adds PAINS/Brenk/aggregator/Ro5 flags and prodrug-active species notes without changing ranking. |
| Upstream docking species resolution | Integrated seam | `vta/nodes/species_resolution.py`, `vta/chem/species.py`, `vta/nodes/docking.py` | Uses curated active-form SMILES when present; otherwise labels known prodrugs/nucleoside analogs as parent-surrogate docking. No active-form structures are fabricated. |
| Catalytic Mg/Mn handling | Integrated seam | `vta/nodes/docking.py` | Docking records consume curated pocket/structure metal coordinates when supplied; otherwise report `[skip] Metals: no curated catalytic Mg/Mn coordinates`. |
| Ensemble docking | Integrated seam, default off | `vta/nodes/docking.py` | Real Vina path can dock explicit receptor ensembles and keep the best conformer; default remains single receptor with a labelled skip. |
| Score provenance | Integrated | `vta/provenance.py` | Adds tool/version/seed/input hash blocks beside score-like fields. |

## What Is Actually Active In A Normal Run

Always active in the graph:

- Classification and routing.
- Structure node with experimental, ESMFold, AlphaFold DB, Boltz-2, or refusal paths.
- ProteinTTT node, but it only acts if the package exists and a low-pLDDT ESMFold fold exists.
- Pocket node with experimental active site, FPocket/P2Rank, or mock fallback.
- Conservation node with MSA-based JSD or neutral fallback.
- Species-resolution node that chooses parent vs. curated active-form docking species.
- Docking node with Vina or deterministic mock fallback.
- Contact-conservation node when real poses and residue maps exist.
- GNINA and Boltzina nodes, but as annotation-only seams.
- LE-led ranking.
- ADMET node, but only annotates when `admet_ai` is installed.
- HTML report generation.

Not active unless explicitly enabled:

- MD validation, unless `vta run --include-md` or `build_app(include_md=True)` is used.
- Production MM-GBSA, because it requires a completed MD run plus a compatible engine.
- FEP/ABFE, unless `vta run --include-fep` or `build_app(include_fep=True)` is used,
  and an external `openfe-abfe` runner is available.

## Integration Pattern Used

Every accuracy feature follows the same pattern:

1. Add a small node or data helper.
2. Detect the dependency at runtime.
3. Write additive fields into `VTAState`.
4. Keep the graph running when the dependency is absent.
5. Record real/skipped/fallback behavior in `audit_trail` and `versions`.
6. Add hermetic tests that monkeypatch heavy tools and network calls.
7. Keep new scientific score terms out of ranking until the validation gate is recalibrated.

This pattern is visible in `admet.py`, `rescore.py`, `boltzina.py`, `proteinttt.py`,
`md_simulate.py`, and `md_analyze.py`.

## Phase 1 Active-Species Status

The committed control-set mappings are mechanistically labelled but do not yet include
curated active-form SMILES:

- Remdesivir parent prodrug -> remdesivir triphosphate (GS-443902).
- Sofosbuvir parent prodrug -> sofosbuvir triphosphate (GS-461203).
- Molnupiravir parent prodrug -> NHC triphosphate.

Because those active-form structures are not committed, the docking node currently keeps
the parent SMILES for these controls and marks the row as `species_source =
parent_surrogate`. If a future curated `active_form_smiles` is added, the same resolver
will prepare that active species before docking. This is intentionally allowed to change
the positive-control gate result; the gate should be reported before/after, not tuned to
preserve the old 3/4 recovery.

Acceptance check on 2026-06-29:

- `python -m pytest tests -q`: 155 passed.
- `scripts/validate_controls.py` with network access for RCSB fetch: AutoDock Vina 1.2.5
  + FPocket, 3/4 positive controls in top 5, verdict PASS.
- Current gate remains a parent-surrogate control benchmark for remdesivir, sofosbuvir,
  and molnupiravir until curated active-form SMILES/structures are committed.

## Phase 2-7 Implementation Status

Phases 2-7 are now implemented as guarded architecture seams:

- Phase 2: property-matched decoy registry, logAUC, scaffold/time split helpers,
  leakage audit, explicit gate-recalibration criteria, and consensus annotations.
- Phase 3: selectivity, resistance-barrier, and prospective prediction-registry
  contracts. Missing panels are labelled skips.
- Phase 4: target-prioritization node with conservation treated as target-level
  context. Compound ranking weights are unchanged.
- Phase 5: ADMET provider interface with antiviral endpoint placeholders and structure
  QC contracts for pLDDT/resolution/Rfree plus MolProbity/PAE/PocketMiner seams.
- Phase 6: Biolink-aligned KG schema, required edge provenance, and span-grounded
  literature extraction guards that reject negated or ungrounded relations.
- Phase 7: host-virus adapter metadata, database license/access registry, responsible-use
  report note, and machine-readable run provenance block.

All new scientific outputs remain annotation-only unless an explicit validation gate
promotes them. No DL, ADMET, selectivity, resistance, MD, or KG term enters `rank.py`.

## Phase 8 Activation Evidence

Phase 8 is activated on TiLV PB1 / 8PSO chain B.

Current evidence artifacts:

- `outputs/phase8/locked_predictions_tilv_pb1.json`: top predictions are locked with
  timestamp and input hashes.
- `outputs/phase8/matched_benchmark_tilv_pb1.json`: real Vina/FPocket benchmark over
  4 controls plus 50 ChEMBL property-matched presumed decoys.
- `outputs/phase8/gate_decision_matched_tilv_pb1.md`: DL/consensus promotion decision.

Matched benchmark result:

- N = 54: 4 actives / 50 decoys.
- BEDROC(alpha=20) = 0.4062.
- EF1% = 13.5.
- logAUC = 0.3883.
- ROC-AUC = 0.785.
- Controls in top 5 = 1/4.

Decision: no ranking change. The result is useful evidence, but it is not
publication-grade because decoys are presumed inactive, the set is still smaller than
30-50 decoys per active, and no leakage-controlled DL comparison has passed.

## Phase 9 Benchmark Maturity

Phase 9 has started with two concrete maturity upgrades:

- Bootstrap confidence intervals for BEDROC, EF1%, logAUC, and ROC-AUC are implemented
  in `vta/eval/metrics.py`.
- Versioned ChEMBL active-set artifacts are generated under `vta/data/phase9_actives/`.

Current status is tracked in `docs/phase9_benchmark_maturity.md`.

Key result: raw HCV NS5B (`CHEMBL5375`) has 52 ChEMBL activity rows, but
stratification shows they must not be pooled into one pocket: 1 NI, 6 NNI, 8 direct
site-unknown NS5B rows, 35 indirect antiviral-assay rows, and 2 off-target HCV drugs.
The practical powered active-site scope is `hcv_ns5b_ni_active_site.json`, derived
from ChEMBL `CHEMBL4296320`, with 43 NI/nucleoside-like active-site actives. TiLV PB1
remains active-limited, and its current CI-bearing benchmark confirms that 4 actives
cannot support a stable enrichment claim.

## Ranking Impact

Currently ranking is still controlled by `vta/nodes/rank.py`:

```text
score = 0.55 * ligand_efficiency
      + 0.35 * binding_affinity
      + 0.10 * conservation
```

The following are not folded into the main score yet:

- ADMET predictions.
- Active-species/prodrug labels.
- Chemistry flags.
- GNINA CNN score or affinity.
- Boltzina score.
- MM-GBSA energy.
- FEP/ABFE free energy.

That is deliberate. Any new ranking term must be validated with
`scripts/validate_controls.py` and benchmarked before it changes lead ordering.

## Remaining Work

High priority:

- Run GNINA and/or Boltzina on a GPU host and compare whether the annotations improve
  positive-control recovery or enrichment.
- Recalibrate `rank.py` only if DL rescoring improves validation.
- Run a real ADMET-AI pass and confirm report columns with production predictions.
- Run an MD smoke test with `include_md=True` and short environment-controlled step counts.

Medium priority:

- Populate real DUD-E or LIT-PCBA benchmark datasets instead of relying on demo decoys.
- Run `scripts/benchmark_enrichment.py` on a larger active/decoy set.
- Verify Boltzina and ProteinTTT APIs against their installed packages before citing them
  as production-grade integrations.

Future/paper-grade:

- Add publication-ready benchmark tables with EF1%, BEDROC, ROC-AUC, and control recovery.
- Run production FEP/ABFE after the MD and DL-rescoring gates are stable.

## Environment Reality

The code is mostly integrated, but real accuracy gains depend on available tools:

- CPU/local friendly: ADMET-AI, P2Rank, conservation, reports, enrichment metrics.
- External binary required: FPocket, Vina, GNINA, MM-GBSA engines.
- GPU or long compute required: GNINA CNN in production, Boltz/Boltzina, ProteinTTT,
  full MD, FEP/ABFE.

The current implementation correctly separates "seam exists" from "real production run
has been executed." A missing binary or package should produce a labelled skip, not an
untrusted fabricated result.

## Suggested Next Verification Order

1. Run the full test suite.
2. Run `scripts/validate_controls.py` with real Vina/FPocket.
3. Run ADMET-AI on the top leads and inspect the HTML report.
4. Run GNINA/Boltzina on a GPU host and save annotations.
5. Re-run positive-control validation with any proposed new score weights.
6. Run short MD smoke test.
7. Run production MD for top leads.
8. Configure an `openfe-abfe` runner and run FEP/ABFE for the final top leads.
