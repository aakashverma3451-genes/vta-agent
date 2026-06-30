# VTA Dependency Matrix

## Core

| Dependency | Purpose | Install Group |
|------------|---------|---------------|
| FastAPI | Future platform API surface | runtime |
| Pydantic | API and service data validation | runtime |
| NetworkX | Graph algorithms and prototype KG operations | runtime |
| LangGraph | Agent/pipeline orchestration | runtime |
| LangChain Core | Message/tool abstractions | runtime |
| Requests | HTTP APIs | runtime |
| jsonschema | TaxonAgent contract validation support | runtime |

## Structural Biology

| Dependency | Purpose | Install Group |
|------------|---------|---------------|
| Biopython | Sequence/structure utilities | runtime |
| RDKit | Molecule parsing, descriptors, conformers | runtime |
| Meeko | Ligand/receptor PDBQT preparation | runtime |
| Gemmi | CIF/PDB conversion | runtime |
| OpenMM | MD simulation | optional `md` |
| MDAnalysis | Trajectory analysis | optional `md` |
| ParmEd | Amber topology export for MM-GBSA | optional `md` |

## Docking

| Dependency | Purpose | Install Group |
|------------|---------|---------------|
| AutoDock Vina | Docking engine | external binary |
| FPocket | Pocket detection | external binary |
| P2Rank | Pocket consensus | external binary |
| GNINA | CNN pose rescoring | external binary |

## Knowledge Graph

| Dependency | Purpose | Install Group |
|------------|---------|---------------|
| Neo4j | Graph database driver | optional `kg` |
| Qdrant | Vector search for literature/evidence | optional `kg` |
| PostgreSQL / psycopg | Relational metadata storage | optional `kg` |

## AI

| Dependency | Purpose | Install Group |
|------------|---------|---------------|
| OpenAI | LLM reasoning providers | optional `ai` |
| Anthropic | LLM reasoning providers | optional `ai` |
| Transformers | Local model and embedding support | optional `ai` |

## ADMET

| Dependency | Purpose | Install Group |
|------------|---------|---------------|
| admet-ai | hERG, bioavailability, solubility predictions | optional `admet` |

## Validation

| Dependency | Purpose | Install Group |
|------------|---------|---------------|
| Pytest | Test runner | dev |
| Coverage | Coverage reports | dev |

## External Tool Policy

External scientific tools remain outside pip dependency resolution. VTA discovers them
through `vta.toolconfig` or the shell `PATH`, and every node must degrade gracefully when
the tool is unavailable.

