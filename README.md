# VTA-Agent — Viral Target Assessment

A LangGraph pipeline that takes a viral genome, classifies it with **TaxonAgent**
(Module 1), and assesses its proteins as drug targets: fold → find pockets → dock →
rank. The expensive biology (pocket detection, docking) is **mocked** in Phase 1 so
the full contract and orchestration are proven before Phase 2 swaps in real tools.

> Full design, diagrams, and the module contract: see [`ARCHITECTURE.md`](ARCHITECTURE.md).

## Pipeline

```
genome.fasta
   │  TaxonAgent.classify()  →  §4.1.1 contract (validated)
   ▼
classify → route_by_confidence ─┬─ proceed/flag → structure → pockets → dock → rank → top-20 leads
                                └─ defer        → return hypotheses to a human
```

## Status (Phase 1)

| Step | Status |
|------|--------|
| Shared state contract (`VTAState`) | ✅ |
| Classify node (real TaxonAgent) | ✅ |
| Confidence router (proceed / flag / defer) | ✅ |
| Structure node — experimental-first + ESMFold (≤400 aa) + chain extraction | ✅ real |
| Ligand library — real ChEMBL compounds (`vta/data/ligands.py`) | ✅ real (Phase 2a) |
| Pocket detection — real FPocket (`vta/nodes/pockets.py`) | ✅ real (Phase 2b) |
| Docking — real AutoDock Vina (`vta/nodes/docking.py`) | ✅ real (Phase 2c) |
| Ranking (composite score) | ✅ real |
| End-to-end on real TiLV genome | ✅ proven |

### Vina setup
`docking_node` auto-detects the Vina engine (`VINA_BIN` or `vina` on PATH) plus the
RDKit/Meeko prep stack; if any is absent it degrades to the mock. On Apple Silicon
there's no pip wheel (Vina needs Boost), so use the official mac x86_64 release binary
— it runs via Rosetta:
```bash
pip install rdkit meeko gemmi          # ligand/receptor PDBQT prep
export VINA_BIN=/path/to/vina          # AutoDock Vina 1.2.5 mac binary
```
Real Vina under emulation is slow, so the real path docks the top pocket per protein.
Validated: remdesivir → PB1 ≈ −6.4 kcal/mol.

### FPocket setup
The node auto-detects fpocket via `FPOCKET_BIN` or `fpocket` on PATH; if absent it
degrades to the mock. On Apple Silicon there's no conda/brew binary — build from
source with `ARCH = MACOSXARM64` in the makefile, then:
```bash
export FPOCKET_BIN=/path/to/fpocket-src/bin/fpocket
```
Note: `conservation` is not produced by fpocket (needs an MSA) — currently a neutral
placeholder. And druggability is low on an isolated apo polymerase subunit, which is
expected; active-site validation against the bound complex is a Phase-2 follow-up.

## What's real vs mocked

- **Real:** classification, contract + schema validation, confidence routing,
  structure folding (experimental PDB / ESMFold), **ligand library (ChEMBL: real
  ChEMBL ids + SMILES, cached to `vta/data/ligands_chembl.json`)**, ranking, audit
  trail.
- **Mocked (Phase 2):** FPocket pocket detection, AutoDock Vina docking *scores*.
  Mocks return real-shaped data over the real ligands, so the swap is invisible
  downstream — Phase 2 only replaces the fabricated ΔG with real Vina output.

## Run

```bash
# uses TaxonAgent's environment (it carries the classifier + deps)
PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python -m pytest tests/ -q
```

## Known Phase-2 work

- Local GPU ESMFold or chunking (public API has a hard 400-residue ceiling, so
  PB1/PB2 can't use it).
- Resolve the PB2 chain mapping (8PSO labels PA-like / putative-PB1 / RdRp — no
  explicit PB2 chain).
- Real FPocket + Vina + a ChEMBL antiviral/RdRp-inhibitor ligand set (the 5 named
  positive controls must rank top-decile).
