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
| Pockets / docking | ⚠️ mocked (real-shaped) — Phase 2 |
| Ranking (composite score) | ✅ real |
| End-to-end on real TiLV genome | ✅ proven |

## What's real vs mocked

- **Real:** classification, contract + schema validation, confidence routing,
  structure folding (experimental PDB / ESMFold), ranking, audit trail.
- **Mocked (Phase 2):** FPocket pocket detection, AutoDock Vina docking, ChEMBL
  ligand library. Mocks return real-shaped data so the swap is invisible downstream.

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
