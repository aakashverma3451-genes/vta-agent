# VTA-Agent Master Document

## 1. Executive Overview

VTA-Agent is an autonomous viral target assessment and antiviral lead-ranking
pipeline. It takes a viral genome FASTA, classifies the virus through TaxonAgent,
selects protein targets, decides whether confidence is high enough to proceed,
builds or fetches target structures, identifies pockets, screens antiviral
compounds, ranks leads, annotates safety/drug-likeness, and writes an HTML report.

Core flow:

```text
genome FASTA -> TaxonAgent contract -> structures -> pockets -> docking
-> ranked antiviral leads -> report + audit trail
```

The pipeline is built to run in imperfect environments. Scientific tools are used
when present, and missing tools produce labelled fallbacks instead of silent fake
results or crashes.

## 2. Project Goal

The goal is to automate early antiviral discovery triage for viral genomes:

- Identify the virus and extract candidate viral protein targets.
- Route low-confidence classifications to human review.
- Use experimental structures first, then prediction/fallback paths.
- Detect catalytic sites or binding pockets.
- Dock a curated antiviral ligand library.
- Rank compounds with a ligand-efficiency-led composite score.
- Annotate top leads with ADMET and optional validation signals.
- Produce a reproducible report with tool versions and audit decisions.

This is a research and prioritization system, not a clinical decision system.

## 3. Technology Stack

- Language: Python 3.10+
- Package config: `pyproject.toml`, setuptools
- Graph orchestration: LangGraph `StateGraph`
- CLI command: `vta run <genome.fasta>`
- Core dependencies: RDKit, Meeko, Gemmi, Requests, jsonschema
- Optional tools: FPocket, P2Rank, AutoDock Vina, GNINA, ADMET-AI, OpenMM,
  MDAnalysis, AmberTools/gmx_MMPBSA, Boltz/Boltzina, ProteinTTT
- Tests: pytest

The installed console script is:

```text
vta = "vta.cli:main"
```

## 4. System Architecture

VTA-Agent is the target-assessment half of a two-module system:

- TaxonAgent answers: "What virus is this genome?"
- VTA-Agent answers: "What target and leads should be prioritized?"

The boundary is a validated TaxonAgent contract returned by
`taxonagent.classify_genome()`. VTA-Agent validates that contract, copies relevant
fields into `VTAState`, and then lets graph nodes progressively enrich that state.

The main architecture pattern is:

```text
VTAState -> node -> VTAState -> node -> VTAState
```

Nodes do not call each other directly. `vta/graph.py` wires the LangGraph edges.

## 5. Runtime Pipeline

Default fast path:

```text
classify -> router -> structure -> proteinttt -> pockets -> conservation
-> dock -> conservation_contacts -> rescore -> boltzina -> rank -> admet -> report
```

Low-confidence path:

```text
classify -> router -> defer -> report
```

Optional MD path:

```text
... -> admet -> md_select -> md_simulate -> md_analyze -> md_rerank -> report
```

The CLI currently calls `build_app()` without MD. MD is available when code calls
`build_app(include_md=True)`.

## 6. Shared State

`vta/state.py` defines `VTAState`, the single state object passed through the graph.
Important fields:

- `genome_fasta`, `run_id`: run input and output label.
- `taxon_result`, `extracted_proteins`: TaxonAgent output and selected targets.
- `classification_confidence`, `classification_basis`, `route`: routing evidence.
- `structures`, `pockets`, `residue_conservation`: target and pocket data.
- `docking_results`, `lead_candidates`: screening and ranking outputs.
- `md_candidates`, `md_results`, `md_analysis`, `md_validated_leads`: optional MD.
- `audit_trail`, `versions`: reproducibility and provenance.

`new_state()` initializes the input, run ID, audit trail, and version map.

## 7. Codebase Map

Root files:

- `README.md`: quickstart, pipeline diagram, example output, and status table.
- `ARCHITECTURE.md`: deeper architecture and TaxonAgent/VTA-Agent contract design.
- `INTEGRATION_PLAN.md`: roadmap for accuracy and scientific integrations.
- `pyproject.toml`: package metadata, dependencies, and CLI entry point.
- `requirements.txt`: environment dependency list.

Main package:

- `vta/cli.py`: CLI, TaxonAgent path detection, run summary.
- `vta/graph.py`: LangGraph node and edge wiring.
- `vta/state.py`: shared state contract.
- `vta/toolconfig.py`: external binary discovery.
- `vta/report.py`: pure HTML rendering and report writing.
- `vta/__main__.py`: `python -m vta` entry point.

Pipeline nodes:

- `vta/nodes/classify.py`: TaxonAgent call and contract validation.
- `vta/nodes/router.py`: `proceed`, `flag`, or `defer` routing.
- `vta/nodes/structure.py`: RCSB, ESMFold, AlphaFold DB, Boltz-2, or refusal.
- `vta/nodes/proteinttt.py`: optional low-confidence fold refinement.
- `vta/nodes/pockets.py`: experimental active site, FPocket/P2Rank, or mock pockets.
- `vta/nodes/conservation.py`: MSA-based Jensen-Shannon conservation.
- `vta/nodes/docking.py`: AutoDock Vina or deterministic mock docking.
- `vta/nodes/conservation_contacts.py`: ligand-contact-weighted conservation.
- `vta/nodes/rescore.py`: optional GNINA CNN annotation.
- `vta/nodes/boltzina.py`: optional Boltzina binding-score annotation.
- `vta/nodes/rank.py`: composite scoring and top-lead selection.
- `vta/nodes/admet.py`: optional ADMET-AI annotations.
- `vta/nodes/report.py`: graph terminal report node.
- `vta/nodes/md_*.py`: optional molecular-dynamics validation phase.
- `vta/nodes/stubs.py`: defer node and legacy structure stub.

Data, evaluation, and scripts:

- `vta/data/ligands.py`: curated antiviral library with ChEMBL/PubChem refresh.
- `vta/data/ligands_chembl.json`: committed ligand cache for offline runs.
- `vta/data/databases.py`: machine-readable database registry.
- `vta/data/msa.py`, `vta/data/uniprot.py`, `vta/data/pubchem.py`: data helpers.
- `vta/eval/metrics.py`: enrichment factor, BEDROC, ROC-AUC, and reports.
- `scripts/validate_controls.py`: known RdRp inhibitor recovery gate.
- `scripts/benchmark_enrichment.py`: enrichment benchmark driver.
- `scripts/build_conservation_msa.py`: offline MSA construction helper.
- `tests/`: hermetic pytest suite for nodes, graph behavior, reports, and metrics.

## 8. Node Behavior

Classification:
`classify_node` calls TaxonAgent, validates the contract, stores extracted proteins,
confidence fields, versions, and an audit entry.

Routing:
`route_by_confidence` applies thresholds: `>=95` proceeds, `85-95` proceeds with a
flag, and `<85` defers. It records whether the confidence number is an amino-acid
identity band rather than a calibrated probability.

Structure:
`structure_node` prefers experimental PDB mappings for known PB1/PA targets. It then
tries ESMFold for proteins up to 400 aa, AlphaFold DB when a UniProt accession is
available, Boltz-2 when installed, or a labelled refusal.

ProteinTTT:
`proteinttt_node` refines only low-confidence ESMFold structures below pLDDT 70. It
does not alter experimental, AlphaFold, or Boltz-2 structures.

Pockets:
`pockets_node` uses known experimental active sites first. Otherwise it runs FPocket,
adds P2Rank consensus if available, or creates deterministic mock pockets.

Conservation:
`conservation_node` computes per-pocket and per-residue Jensen-Shannon conservation
from a homolog MSA. Without an MSA, it keeps a neutral `0.5` placeholder and logs it.

Docking:
`docking_node` runs real Vina when Vina plus RDKit/Meeko preparation are available.
It caches prepared ligands and writes poses under `structures/`. Receptor prep is
Meeko-first with an OpenBabel fallback (`_prep_receptor`): OpenBabel produces a rigid
receptor PDBQT for structures Meeko 0.7.1 declines (e.g. every SARS-CoV-2 Mpro chain),
leaving the working Meeko path untouched. Without real tools, it creates deterministic
mock docking records over the ligand library.

Contact conservation:
`conservation_contacts_node` uses real pose contacts to make conservation
ligand-specific. Mock records or missing MSAs keep pocket-level conservation.

Annotation nodes:
`rescore_node`, `boltzina_node`, and `admet_node` add extra fields when their tools are
installed. These are annotation-only and do not change ranking today.

Ranking:
`rank_node` scores docking records with ligand efficiency, binding affinity, and
conservation:

```text
score = 0.55 * LE + 0.35 * affinity + 0.10 * conservation
```

Ligand efficiency leads the score to reduce raw docking size bias.

Reporting:
`report_node` writes `outputs/{run_id}_report.html` through `vta.report.write_report`.
The report includes classification, structures, ranked leads, optional MD results,
provenance, versions, and the full audit trail.

## 9. External Tool Discovery

`vta/toolconfig.py` resolves tools in this order:

- Explicit environment variable, such as `VINA_BIN` or `FPOCKET_BIN`.
- Binary on `PATH`.
- Known local build locations under the repo or workspace.

This discovery layer is why the graph can run in CI with fallbacks and on a fully
configured machine with real tools.

## 10. Outputs

- `outputs/`: generated HTML reports and benchmark outputs.
- `outputs/{run_id}/md/`: optional MD trajectories, logs, and topology files.
- `structures/`: fetched/generated PDBs, prepared receptors, ligand PDBQT cache,
  Vina poses, and pocket-tool outputs.
- `validation.log`: validation run log.

Many files under `outputs/` and `structures/` are generated artifacts from previous
runs, not core source files.

## 11. Running The Project

Install:

```bash
pip install -e .
```

Run:

```bash
vta run path/to/genome.fasta
```

Test:

```bash
python -m pytest tests -q
```

TaxonAgent is expected to be installed as a package, either from the sibling checkout
with `pip install -e ../taxonagent` or from a package index.

## 12. Current Status

- Scientific validation: **Phases 0–9 complete in software** (see Section 15). Phase 9
  meets its acceptance bar: ≥20 actives, 3 targets, every benchmarkable metric with a 95%
  CI, an executed parent-vs-active-form comparison, and a frozen benchmark artifact.
- Real paths: TaxonAgent seam, routing, graph, report generation, ligand library,
  ranking, conservation logic, enrichment metrics, and several tool integration seams.
- Tool-dependent real paths: RCSB/ESMFold/AlphaFold/Boltz-2, FPocket/P2Rank, Vina,
  GNINA, ADMET-AI, OpenMM/MDAnalysis/MM-GBSA.
- Fallback-capable paths: structures, pockets, docking, conservation, rescoring,
  ADMET, and MD.
- Annotation-only today: GNINA, Boltzina, ADMET, and MM-GBSA.
- Optional slow path: MD validation.

## 13. Validation Strategy

The project guards against raw docking bias with positive controls and enrichment
metrics:

- Known RdRp inhibitors are flagged in the ligand library.
- `scripts/validate_controls.py` checks control recovery.
- `vta/eval/metrics.py` provides EF, BEDROC, and ROC-AUC.
- Ranking emphasizes ligand efficiency because raw docking affinity tends to reward
  larger molecules.

## 14. Design Principles

- Keep TaxonAgent and VTA-Agent separated by one validated contract.
- Pass all graph data through one typed state object.
- Defer when classification evidence is weak.
- Prefer real scientific tools, but label missing-tool fallbacks honestly.
- Keep heavy model outputs as annotations until validation justifies ranking changes.
- Preserve provenance through versions and the audit trail.

## 15. Validation Results (Phases 8–9, current)

The validation is deliberately built to deflate its own inflation. Three targets, every
benchmarkable metric as median + 95% bootstrap CI. Frozen artifact:
`outputs/phase9/locked_benchmark.json` (hash-stamped, immutable). Full narrative:
`docs/phase9_benchmark_maturity.md`; multi-target table:
`outputs/phase9/multitarget_ci_table.md`.

- **TiLV PB1** (8PSO:B) — underpowered demonstration. 4 actives + 50 property-matched
  decoys. A 9-ligand gate showed 3/4 controls in top-5; with matched decoys this collapsed
  to 1/4 and bootstrap CIs are ≈ [0,1]. Reported as a demonstration, not a claim.
- **HCV NS5B NI** (43 nucleotide/triphosphate actives) — un-benchmarkable by matched decoys
  (the nucleotide meta-finding). 24/43 actives recover **zero** scaffold-distinct
  property-matched decoys because triphosphates have no inactive property twins.
  Matched-decoy validation systematically cannot assess nucleotide-analog antivirals;
  pipelines that appear to are usually docking the parent prodrug.
- **SARS-CoV-2 Mpro (non-covalent)** — the powered target. Real AutoDock Vina, 50
  non-covalent COVID-Moonshot actives vs 50 experimentally measured Moonshot inactives,
  7L11:A. BEDROC(α=20) 0.68 [0.36, 0.89], logAUC 0.20 [0.14, 0.30], ROC-AUC 0.58
  [0.47, 0.69], EF1% 2.0 (ceiling-limited at 50/50). Benchmark **design** is
  publication-grade (real measured inactives, binding-mode stratified, informative CIs —
  unlike TiLV's [0,1]); the rigid-Vina docking **signal** is modest (ROC-AUC CI crosses
  0.5) and reported as-is.
- **Phase 9D (parent vs active-form)** — 17 parent prodrugs vs 26 curated triphosphate
  active forms docked into NS5B 2XI3:A with experimental Mg²⁺: median ΔG difference only
  −0.10 kcal/mol (sofosbuvir parent even out-scored its own triphosphate). Poor nucleotide
  recovery is a genuine Vina scoring limitation, not a prodrug-input artifact.
- **Gate** — DO NOT PROMOTE. No DL/consensus term beats the Vina baseline at non-overlapping
  CIs (no DL rescore wired); consensus/DL stays annotation-only, ranking weights unchanged.

Supporting infrastructure: covalent-warhead binding-mode stratification
(`vta/chem/warheads.py`); Meeko-first receptor prep with OpenBabel fallback
(`vta/nodes/docking._prep_receptor`); bootstrap + multi-draw CIs in `vta/eval/metrics.py`.
Phase-9 drivers: `scripts/fetch_mpro_dataset.py`, `scripts/stratify_mpro_binding_mode.py`,
`scripts/run_mpro_benchmark.py`, `scripts/run_phase9d_active_form.py`,
`scripts/build_multitarget_summary.py`, `scripts/freeze_benchmark.py`.

### Phase 11 — validity hardening (the peer-review layer)

The most consequential result: **on the powered Mpro benchmark, AutoDock Vina does NOT beat a
trivial 2D-similarity baseline** (2D-sim BEDROC 0.92 vs Vina 0.68; paired-bootstrap Δ 95% CI
[−0.57, +0.06]) — so structure-based enrichment is *not demonstrated* (Wallach & Heifets 2018),
reported honestly. Diagnostics show this is a ranking/benchmark limit, not a protocol failure:
Vina redocks the native Mpro ligand to 1.65 Å (pose-reliable), and the 9D parent-vs-active-form
ΔG gap is only −0.10 kcal/mol. Method upgrades: trivial baselines (`vta/eval/baselines.py`),
paired-bootstrap significance replacing all CI-overlap logic (`vta/eval/significance.py`), a
gated **LE→ΔG ranking demotion** (`vta/nodes/rank.py` now ΔG-primary), and a GNINA rescorer
that no-ops cleanly to DO-NOT-PROMOTE when absent (`vta/nodes/rescore.py`). The nucleotide
meta-finding was reframed from "impossible" to recipe-specific.

**Extended docks (Phase 11 WI-4/WI-5), closing the last reviewer objections:**

- **Job A — Mpro at a non-degenerate 31:1 ratio** (24 non-covalent actives vs 740 real
  measured Moonshot inactives; EF1% ceiling lifted 2.0 → 31.8;
  `outputs/phase11/mpro_30to1_benchmark.json`). The "Vina < 2D-similarity" finding does not
  just survive the honest ratio — it **strengthens into statistical significance**:
  2D-similarity BEDROC 0.40 [0.24, 0.56] / ROC-AUC 0.82 [0.72, 0.90] vs Vina BEDROC 0.15
  [0.06, 0.28] / ROC-AUC 0.63 [0.51, 0.75]; paired Δ(Vina−2Dsim) BEDROC −0.24, **95% CI
  [−0.44, −0.04]** (P=0.01) and ROC-AUC −0.18, **95% CI [−0.31, −0.07]** (P=0.001) — both
  strictly < 0, so **2D-memorization significantly beats docking**. Vina's EF1% is **0.0**
  (its top ~1% of ranks contains no actives). At the degenerate 1:1 ratio the loss was
  borderline (CI touched 0); at 31:1 it is unambiguous. The docking loop was made
  resumable + hang-capped (180 s/dock) after a few large/flexible inactives stalled Vina.
- **Job B — HCV NS5B NI with property-UNMATCHED decoys** (`scripts/phase11_hcv_unmatched.py`,
  640 ChEMBL drug-like phosphate-free decoys; DeepCoy absent): tests whether the nucleotide
  benchmark runs under a different decoy recipe, and whether any "enrichment" is binding
  signal or mere charge/size discrimination. **[result pending — dock in progress]**

A slide-ready narrative lives in `docs/VTA_AGENT_DECK.md`.
