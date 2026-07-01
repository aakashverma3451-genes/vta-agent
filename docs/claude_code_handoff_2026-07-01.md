# Claude Code Handoff — VTA-Agent Phase 10 (2026-07-01)

Supersedes `docs/claude_code_handoff_2026-06-30.md`. Branch: `session-a`.

## Current Position

**Phase 10 (ensemble docking) complete.** Phase 9 + 10 form a complete
mechanistic-correctness block: (a) single-structure powered benchmark with CIs,
(b) parent-vs-active-form comparison ruling out wrong-species artifact (Phase 9D),
(c) 3-conformer receptor ensemble quantifying conformational stability (Phase 10).

## Phase 10 Results

Ensemble: 7L11:A (primary holo) + 6Y2E_apo (6Y2E:A), 7K3T_holo (7K3T:A).
N compounds: 100 (same Phase 9 set).
New docking runs: 2 conformers × 100 compounds = 200 calls.

| Metric | Phase 9 (7L11:A) | Phase 10 (ensemble) | Delta | CI overlap? |
|--------|-----------------|---------------------|-------|-------------|
| BEDROC(α=20) | 0.6767 [0.357, 0.893] | 0.7848 [0.461, 0.932] | +0.1081 | True |
| ROC-AUC | 0.5796 [0.467, 0.690] | 0.4686 [0.347, 0.578] | -0.1110 | True |

Top-10 Phase 9 actives with improved ensemble rank: 3/10.

**Interpretation:** Ensemble (3 conformers) BEDROC 0.7848 vs Phase 9 single-structure 0.6767 (delta +0.1081). CI bands overlap — ensemble change is within noise; Phase 9 single-structure signal is stable w.r.t. conformational sampling. Top-10 Phase 9 actives: 3/10 improved rank in ensemble. Gate unchanged: DO NOT PROMOTE without leakage-controlled DL improvement at non-overlapping CIs.

## Gate (unchanged)

DO NOT PROMOTE — no DL rescore wired. Ensemble conformer analysis does not
itself license a ranking change. See `outputs/phase9/gate_decision_mpro.md`.

## Test Status

`PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python -m pytest tests -q`
→ **215 passed.**

## Non-Negotiable Constraints (carry forward)
- No fabricated SMILES / decoys / metal coordinates / assay data.
- Never pool binding modes (covalent vs non-covalent; NI vs NNI).
- Every enrichment metric reported as median + 95% CI, never a bare point.
- No ranking change outside the gate; weights unchanged.
- Do not inflate the modest Mpro docking signal.

## Commit Status
All Phase 10 work committed on session-a.

## Next Work — Phase 12 (API / Deployment)
- FastAPI endpoint + async queue (Celery/Redis or asyncio) for pipeline runs
- Persistence layer (SQLite or Postgres) for run history + result retrieval
- Containerization (Docker + docker-compose for Vina + pipeline deps)
- Phase 11 (wet-lab) is evidence-gated parallel track — lock predictions,
  ingest assay results when available.
- Phase 13 (RO-Crate/PROV-O, Zenodo, manuscript) after Phase 12.

## Commands
```bash
# tests
PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python -m pytest tests -q
# Phase 10 ensemble (rerun, ~2-3h for 200 new dockings)
PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python scripts/run_phase10_ensemble.py
```

## Recommended Next Claude Code Prompt
```text
Continue VTA-Agent from docs/claude_code_handoff_2026-07-01.md. Start Phase 12:
FastAPI endpoint + async job queue for pipeline runs, SQLite persistence,
Dockerfile for Vina + pipeline deps. Keep all tests green. No fabrication, no
ranking change outside the gate.
```
