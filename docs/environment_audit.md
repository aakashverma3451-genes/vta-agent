# VTA Environment Audit

## Current State

VTA-Agent previously relied on the sibling TaxonAgent virtual environment for test and
script execution. The common local command was:

```bash
PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python -m pytest tests -q
```

That made VTA operational, but it coupled VTA to:

- The sibling `../taxonagent` checkout.
- The sibling `../taxonagent/venv` interpreter.
- A manual `PYTHONPATH` value.
- The internal TaxonAgent source layout.

Runtime coupling has now been removed from `vta.cli`: VTA imports `taxonagent` as a
normal installed package. Development can still use `pip install -e ../taxonagent`, but
VTA no longer prepends `../taxonagent/src` at runtime.

## Current Python Environments

Observed in this workspace:

- System `python3`: available, but does not have pytest installed.
- `../taxonagent/venv/bin/python`: available and currently used for validation.
- No project-local `.venv` or conda environment was present in this repository.

Target environment:

- Conda/mamba environment named `vta`.
- Python 3.11.
- VTA installed with `pip install -e .`.
- TaxonAgent installed with `pip install -e ../taxonagent` or `pip install taxonagent`.

## Dependency Graph

```text
VTA-Agent
├── taxonagent
│   └── classify_genome / validate_contract
├── orchestration
│   ├── langgraph
│   └── langchain-core
├── chemistry + structure
│   ├── rdkit
│   ├── meeko
│   ├── gemmi
│   ├── biopython
│   ├── openmm
│   └── mdanalysis
├── external binaries
│   ├── fpocket
│   ├── vina
│   ├── gnina
│   └── openfe-abfe
├── optional models
│   ├── admet-ai
│   ├── boltzina
│   └── proteinttt
├── platform services
│   ├── fastapi
│   ├── pydantic
│   ├── neo4j
│   ├── qdrant-client
│   └── psycopg
└── AI clients
    ├── openai
    ├── anthropic
    └── transformers
```

## Shared Dependencies

Likely shared with TaxonAgent:

- `requests`
- `jsonschema`
- pytest/dev tooling
- Possibly sequence-biology packages, depending on TaxonAgent internals

VTA-specific:

- RDKit, Meeko, Gemmi
- LangGraph orchestration
- Vina/FPocket/GNINA seams
- OpenMM/MDAnalysis/ParmEd
- ADMET, knowledge graph, vector DB, and AI client dependencies

## Potential Conflicts

- RDKit is more reliable from conda-forge than pip on many platforms.
- OpenMM, OpenMMForceFields, and OpenFF Toolkit can be sensitive to Python and conda
  channel versions.
- Torch-backed packages such as `admet-ai`, `transformers`, Boltzina, and ProteinTTT
  may pull large or platform-specific wheels.
- External binaries are not Python dependencies and must remain autodetected.
- TaxonAgent must expose `classify_genome` and `validate_contract` at package import
  time; older installed versions may be incompatible.

## Recommended Separation Plan

1. Create the standalone conda environment:

```bash
mamba env create -f environment.yml
mamba activate vta
```

2. Install TaxonAgent as a package:

```bash
pip install -e ../taxonagent
```

or:

```bash
pip install taxonagent
```

3. Install VTA:

```bash
pip install -e ".[dev]"
```

4. Keep external scientific binaries outside Python dependency resolution:

- `vina`
- `fpocket`
- `gnina`
- `openfe-abfe`

5. Run validation without sibling interpreter coupling:

```bash
python -m pytest tests -q
```

## Import Path Target

Target state:

```text
TaxonAgent
  -> installed package
  -> imported by vta.nodes.classify
  -> consumed by VTA graph
```

VTA should not require:

- `../taxonagent/venv`
- `../taxonagent/src`
- manual `sys.path` mutation
- runtime `PYTHONPATH`

