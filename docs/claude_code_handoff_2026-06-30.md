# Claude Code Handoff — VTA-Agent Scientific Validation (2026-06-30)

Supersedes `docs/claude_code_handoff_2026-06-29.md`. Workspace:
`/Users/dna/Desktop/vtaagent/vta-agent`. Branch: `session-a`.

## Current Position

Phases 0–9 are complete in software. **Phase 9 meets its acceptance bar** (≥20 actives,
≥2 targets, every metric with a 95% CI, an executed parent-vs-active-form comparison, and a
frozen benchmark artifact). The honest core of the project is intact and now stronger: the
pipeline is built to deflate its own inflation, and two of the three target results are
negative/limited — framed as the contribution, not hidden.

## Headline Results (three targets, all CI-bearing where benchmarkable)

Frozen in `outputs/phase9/locked_benchmark.json` (hash-stamped; immutable — new dated
artifact on change). Multi-target table: `outputs/phase9/multitarget_ci_table.md`.

- **TiLV PB1** (8PSO:B, ChEMBL property-matched presumed decoys): underpowered demonstration,
  4 actives. BEDROC 0.41, ROC-AUC ~0.79, but bootstrap CIs ≈ [0,1]. Not a claim.
- **HCV NS5B NI / triphosphate** (43 actives): **un-benchmarkable by matched decoys** — the
  nucleotide meta-finding. 24/43 actives recover **zero** scaffold-distinct property-matched
  decoys (3.09/avg); triphosphates have no inactive property twins.
- **SARS-CoV-2 Mpro (non-covalent)** — the powered target. Real AutoDock Vina, 50 non-covalent
  Moonshot actives vs 50 **experimentally measured** Moonshot inactives, 7L11:A active site:
  - BEDROC(α=20) 0.68, 95% CI **[0.357, 0.893]**
  - logAUC 0.20, CI **[0.137, 0.301]**
  - ROC-AUC 0.58, CI **[0.467, 0.690]**
  - EF1% 2.0 (ceiling-limited at 50/50)
  - **Two honest layers:** benchmark *design* is publication-grade (real measured inactives,
    binding-mode stratified, real Vina, bootstrap + multi-draw CIs, informative intervals —
    unlike TiLV's [0,1]); docking *signal* is **modest** (ROC-AUC CI crosses 0.5). Reported
    as-is, not inflated.

### Phase 9D — parent vs active-form (executed, the tie-breaker)
Docked 17 parent prodrugs vs 26 ChEMBL-curated triphosphate active forms into the real HCV
NS5B catalytic site (**2XI3:A, experimental catalytic Mg²⁺ retained**, box on the GDD
asp-220/318/319 centroid). `outputs/phase9/phase9d_active_form_comparison.json`.
- Median ΔG: triphosphate −7.57 vs parent −7.47 → **difference −0.10 kcal/mol** (within Vina
  error). Matched pair: sofosbuvir parent (−7.90) scored *stronger* than its triphosphate
  GS-461203 (−7.35).
- **Conclusion:** poor nucleotide recovery is **not** primarily a "wrong species docked"
  artifact — rigid Vina does not reward the triphosphate's catalytic-metal coordination. It
  is a genuine scoring limitation, which mechanistically explains the un-benchmarkability.
- TiLV active-form re-dock stays a **labelled skip** (no curated triphosphate SMILES; not
  fabricated).

### Gate
`outputs/phase9/gate_decision_mpro.md`: **DO NOT PROMOTE.** Even a powered benchmark does not
license a ranking change without a leakage-controlled DL improvement at non-overlapping CIs,
and no DL rescore is wired (`consensus_node` `cnn_affinity`/`boltzina_score` unpopulated).

## Infrastructure added since 2026-06-29

- **OpenBabel receptor-prep fallback** in `vta/nodes/docking._prep_receptor`: Meeko-first
  (TiLV path untouched), then a rigid OpenBabel PDBQT when Meeko declines. This unblocked
  Mpro, where Meeko 0.7.1 fails deterministically (`update_H_positions`) on every Mpro chain.
  Dep added: `openbabel-wheel` in `requirements.txt`.
- **Covalent-warhead classifier** `vta/chem/warheads.py` (Cys145 electrophiles) — used to
  hold covalent Mpro inhibitors out of the non-covalent benchmark (binding-mode discipline,
  same rule that separated HCV NI from NNI).
- **Mpro dataset pipeline**: `scripts/fetch_mpro_dataset.py` (Moonshot CSV primary + ChEMBL
  CHEMBL4523582 top-up), `scripts/stratify_mpro_binding_mode.py`. Data in `vta/data/mpro/`.
- **Benchmark + freeze**: `scripts/run_mpro_benchmark.py` (Vina-affinity enrichment, bootstrap
  + multi-draw CIs), `scripts/run_phase9d_active_form.py`, `scripts/freeze_benchmark.py`,
  `scripts/build_multitarget_summary.py`.
- Mpro wired into `EXPERIMENTAL_PDB`/`EXPERIMENTAL_ACTIVE_SITE` (7L11:A, XF1 centroid).

## Test Status
`PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python -m pytest tests -q` → **200 passed.**

## Non-Negotiable Constraints (carry forward)
- No fabricated SMILES / decoys / metal coordinates / assay data — missing → labelled skip
  with provenance.
- Never pool binding modes (covalent vs non-covalent; NI vs NNI; NI-site vs allosteric).
- Every enrichment metric reported as median + 95% CI, never a bare point.
- No ranking change outside the gate; weights unchanged.
- Do not inflate the modest Mpro docking signal; do not claim TiLV enrichment.
- Leave the `.claude/COLLAB.md` deletion alone.

## Commit Status (IMPORTANT)
- Committed: `bce4182` — Phase 9C scaffolding + the prior research-accuracy foundation bundle.
- **Uncommitted (14 paths):** the OpenBabel fallback (`vta/nodes/docking.py`,
  `requirements.txt`), the real Mpro benchmark + updated table/gate/doc, and all of Phase
  9D/9E (`run_phase9d_active_form.py`, `freeze_benchmark.py`, `test_phase9d_active_form.py`,
  `phase9d_active_form_comparison.*`, `locked_benchmark.json`). Commit before continuing.

## Next Work — Phase 10 (mechanistic correctness, executed)
- Ensemble docking on top leads (apo/holo + MD/AlphaFold/Boltz-2 conformers); quantify whether
  ensembles change ranking vs single-structure.
- Finalize protonation/tautomer/metal handling as the default for the validation targets.
- Done when: the Phase 9 benchmark is regenerated on ensemble + correct-protomer inputs with
  the change quantified.
- Then Phase 12 (FastAPI/queue/persistence/containerization) and Phase 13 (RO-Crate/PROV-O,
  Zenodo splits+decoys+actives, manuscript from generated artifacts). Phase 11 (wet-lab) is
  the parallel evidence-gated track (lock + ingest predictions; VTA does not design assays).

## Commands
```bash
# tests
PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python -m pytest tests -q
# Mpro powered benchmark (real Vina; ~50-90 min for 50+50)
PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python scripts/run_mpro_benchmark.py
# parent vs active-form (HCV NS5B, ~20-30 min)
PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python scripts/run_phase9d_active_form.py
# refresh multi-target table + gate, then re-freeze
PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python scripts/build_multitarget_summary.py
PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python scripts/freeze_benchmark.py
```

## Recommended Next Claude Code Prompt
```text
Continue VTA-Agent from docs/claude_code_handoff_2026-06-30.md. First commit the uncommitted
Phase 9C/9D/9E + OpenBabel-fallback work on session-a. Then start Phase 10: ensemble docking
on top Mpro leads and finalize protomer/metal handling, regenerate the Phase 9 Mpro benchmark
on ensemble inputs, and quantify the ranking change vs single-structure. Keep CIs on every
metric, no fabrication, no ranking change outside the gate, tests green.
```
