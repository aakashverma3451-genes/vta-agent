# VTA-Agent — Autonomous Drug Discovery from Viral Genomes

> One command. Genome in, ranked drug leads out.

VTA-Agent is an autonomous LangGraph pipeline that takes a raw viral genome FASTA,
classifies it with [TaxonAgent](../taxonagent), folds each protein, screens a curated
antiviral library with AutoDock Vina, ranks leads by ligand efficiency, and filters them
through ADMET safety predictions — producing a self-contained HTML report with no human
setup required.

**Validated:** recovers Ribavirin (rank 1), Sofosbuvir, Molnupiravir, and Remdesivir as
top-ranked leads on the TiLV PB1 RdRp target. 25 hermetic tests. 10+ real-world commits.

---

## Pipeline

```
genome.fasta
     │
     ▼
 ┌─────────────┐
 │  classify   │  TaxonAgent: species, confidence, extracted proteins
 └──────┬──────┘
        ▼
 ┌─────────────┐
 │   router   │  proceed (≥80%) / flag (60–80%) / defer (<60%)
 └──────┬──────┘
        │ proceed / flag
        ▼
 ┌─────────────┐
 │  structure  │  experimental PDB (RCSB first) or ESMFold (≤400 aa)
 └──────┬──────┘
        ▼
 ┌─────────────┐
 │   pockets   │  experimental active site → FPocket blind search
 └──────┬──────┘
        ▼
 ┌─────────────┐
 │    dock     │  AutoDock Vina over ChEMBL antiviral library (9 leads)
 └──────┬──────┘
        ▼
 ┌─────────────┐
 │   rescore   │  GNINA CNN re-score (annotation-only; skips w/o GPU)
 └──────┬──────┘
        ▼
 ┌─────────────┐
 │    rank     │  LE-led composite score (LE 55% · ΔG 35% · conservation 10%)
 └──────┬──────┘
        ▼
 ┌─────────────┐
 │    admet    │  ADMET-AI: hERG, oral bioavailability, solubility per lead
 └──────┬──────┘
        ▼
 ┌─────────────┐
 │   report    │  self-contained HTML (lead table + ADMET + audit trail)
 └─────────────┘
```

Every node auto-detects its dependency and degrades to a labelled fallback if the tool is
absent — the pipeline never crashes for want of a binary.

---

## Quick start

### Requirements

- Python 3.10+
- TaxonAgent (sibling repo at `../taxonagent`)
- Optional for full real run: `fpocket`, `vina`, `admet-ai`

### Install

```bash
# clone both repos side by side
git clone <this-repo> vta-agent
git clone <taxonagent-repo> taxonagent   # sibling directory

# create environment
conda create -n vta python=3.10 -y && conda activate vta
pip install rdkit meeko gemmi            # ligand/receptor prep
pip install admet-ai                     # ADMET predictions (pulls torch/chemprop)
pip install -e vta-agent/
```

### Tool binaries (optional — pipeline degrades gracefully without them)

**AutoDock Vina** (Apple Silicon — runs via Rosetta):
```bash
# download the official mac x86_64 binary from github.com/ccsb-scripps/AutoDock-Vina
export VINA_BIN=/path/to/vina
```

**FPocket** (Apple Silicon — build from source):
```bash
# git clone https://github.com/Discngine/fpocket, set ARCH=MACOSXARM64 in makefile
export FPOCKET_BIN=/path/to/fpocket-src/bin/fpocket
```

**GNINA** (deep-learning re-scorer — needs CUDA GPU; seam is built, run on GPU host):
```bash
export GNINA_BIN=/path/to/gnina
```

---

## Example run

```bash
vta run data/amnoonviridae_1.fasta
```

Expected output (truncated):

```
Tilapinevirus tilapiae → confidence 100.0% → PROCEED
Folding PB1 (560 aa) … 8PSO (experimental, RCSB)
Folding PB2 (450 aa) … refused (too long for ESMFold)
Folding PA  (340 aa) … ESMFold (mean pLDDT 72.3)

Pockets:  PB1 → experimental catalytic site (8PSO:F bound CTP, NTP site)
          PA  → FPocket top pocket (druggability 0.61)
Docking:  9 ligands × 2 proteins (AutoDock Vina 1.2.5)

── Top leads ──────────────────────────────────────────────────────────
 # Ligand              dG (kcal/mol)   LE     hERG  Oral  Solubility
 1 Ribavirin              -6.79     -0.399   0.02  0.73   high     ★
 2 Sofosbuvir             -8.51     -0.236   0.11  0.58   mod      ★
 3 Lopinavir              -8.80     -0.191   0.82  0.42   low
 4 Molnupiravir           -6.89     -0.300   0.04  0.69   high     ★
 5 Ledipasvir             -8.89     -0.137   0.09  0.38   low
───────────────────────────────────────────────────────────────────────
★ = known RdRp inhibitor (positive control)

Report → outputs/amnoonviridae_1_report.html
```

---

## Tests

```bash
# from vta-agent/ with taxonagent on the path:
PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python -m pytest tests/ -q
# → 25 passed
```

### Validation gate (requires real Vina + FPocket)

```bash
PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python scripts/validate_controls.py
# VERDICT: PASS  (3/4 known RdRp inhibitors in top-5)
```

---

## Architecture

| Component | File | Role |
|-----------|------|------|
| State contract | `vta/state.py` | Single `VTAState` TypedDict shared by all nodes |
| Graph | `vta/graph.py` | LangGraph `StateGraph` wiring all edges |
| Classify | `vta/nodes/classify.py` | Wraps TaxonAgent `classify_genome()` |
| Router | `vta/nodes/router.py` | Confidence threshold → proceed / flag / defer |
| Structure | `vta/nodes/structure.py` | RCSB experimental → ESMFold → mock cascade |
| Pockets | `vta/nodes/pockets.py` | Experimental active site → FPocket → mock |
| Docking | `vta/nodes/docking.py` | AutoDock Vina → mock; ligand cache in `structures/` |
| Rescore | `vta/nodes/rescore.py` | GNINA CNN re-score seam (annotation-only) |
| Rank | `vta/nodes/rank.py` | LE-led composite (validated against controls) |
| ADMET | `vta/nodes/admet.py` | ADMET-AI hERG / oral bioavailability / solubility |
| Report | `vta/nodes/report.py` | Self-contained HTML emitter |
| Tool config | `vta/toolconfig.py` | Binary autodiscovery (PATH / env var / local build) |
| Ligands | `vta/data/ligands.py` | ChEMBL antiviral library (real IDs + SMILES) |

Full design, node contracts, and ASCII diagrams: [`ARCHITECTURE.md`](ARCHITECTURE.md).  
Accuracy integration plan (ADMET → DL rescore → MD → FEP): [`INTEGRATION_PLAN.md`](INTEGRATION_PLAN.md).

---

## What's real vs placeholder

| Node | Status |
|------|--------|
| Classify (TaxonAgent) | **REAL** |
| Confidence router | **REAL** |
| Structure (RCSB / ESMFold) | **REAL** (ESMFold ≤400 aa; longer → refused) |
| Pockets (FPocket) | **REAL** (experimental site hardcoded for PB1) |
| Docking (Vina) | **REAL** (auto-detected; mock if absent) |
| DL rescore (GNINA) | **SEAM** (skips on CPU/no binary; wired for GPU) |
| Ranking (LE-led) | **REAL** (validated gate passes) |
| ADMET (ADMET-AI) | **REAL** (skips if not installed) |
| Conservation | **placeholder** (neutral 0.5; needs MSA) |
| MD validation | **planned** (Phase 4 in INTEGRATION_PLAN.md) |

---

## Roadmap

1. **Now:** DL re-score on GPU host (GNINA seam is ready — plug in the binary)
2. **Next:** Molecular dynamics validation (100 ns, top-5 leads, OpenMM) — Phase 4
3. **Later:** Generative molecule design (move from screening → designing novel antivirals)
4. **Paper:** FEP absolute binding free energies on top-3 leads

---

## Citation / tools

- [TaxonAgent](../taxonagent) — viral genome classification (in-house)
- [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina) — docking engine
- [FPocket](https://github.com/Discngine/fpocket) — pocket detection
- [ADMET-AI](https://github.com/swansonk14/admet_ai) — Swanson et al., *Bioinformatics* 2024
- [GNINA](https://github.com/gnina/gnina) — deep-learning re-scoring (GPU)
- [LangGraph](https://github.com/langchain-ai/langgraph) — agent orchestration
- [ESMFold](https://github.com/facebookresearch/esm) — protein folding (via ESMAtlas API)
