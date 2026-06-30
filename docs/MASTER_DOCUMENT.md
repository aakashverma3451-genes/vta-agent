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
It caches prepared ligands and writes poses under `structures/`. Without real tools,
it creates deterministic mock docking records over the ligand library.

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
