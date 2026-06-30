# VTA Environment Migration Plan

## Phase 1: Environment Separation

- Create independent `vta` conda/mamba environment from `environment.yml`.
- Install TaxonAgent as a package with `pip install -e ../taxonagent` or
  `pip install taxonagent`.
- Install VTA with `pip install -e ".[dev]"`.
- Run tests using the VTA environment's `python`, not `../taxonagent/venv/bin/python`.

## Phase 2: Package Cleanup

- Remove runtime `PYTHONPATH` assumptions.
- Keep TaxonAgent access limited to `vta.nodes.classify`.
- Move future package code from `src/vta` into the production package layout when the
  project formally migrates to `src/`.
- Keep docking, ranking, and positive-control behavior unchanged during migration.

## Phase 3: Knowledge Graph Integration

- Use `src/vta/kg/schema.py` as the initial schema contract.
- Add Neo4j-backed repository implementations behind interfaces.
- Add Qdrant-backed literature/evidence vector search.
- Keep graph ingestion separate from the default docking/ranking pipeline until
  validation criteria are defined.

## Phase 4: Scientific Agents

- Use `src/vta/future/agents.py` contracts for:
  - Drug Knowledge Graph service
  - Literature Agent
  - Target Prioritization Agent
  - Drug Repurposing Agent
  - Experiment Design Agent
- Require evidence records and confidence fields for every agent output.
- Do not allow agents to change scoring defaults without benchmark validation.

## Phase 5: Production Deployment

- Add FastAPI service layer.
- Add persistent run storage.
- Add queue-backed execution for Vina, GNINA, MD, and FEP jobs.
- Add deployment health checks for external binaries and model packages.
- Add reproducibility metadata to every report.

## Migration Risks

- TaxonAgent package version mismatch.
- RDKit/OpenMM dependency solver conflicts.
- GPU-only packages failing on CPU hosts.
- External binary path drift.
- Accidentally changing validated docking/ranking behavior during packaging cleanup.

## Acceptance Criteria

- `python -m pytest tests -q` passes inside the standalone `vta` environment.
- `vta run <genome.fasta>` works without `PYTHONPATH`.
- `scripts/validate_controls.py` still passes with real Vina/FPocket.
- No ranking or docking behavior changes are introduced by environment migration.

