# VTA-Agent — Multi-Session Collaboration Board

More than one Claude Code session may run against this repo at once. Coordinate here + git.
Rules: read this at task start; claim files under ACTIVE WORK before editing (add `LOCKED`);
move to DONE/HANDOFF when finished; commit small on your branch.

---

## ACTIVE WORK

- **session-a (fable)** — 2026-06-30 — **Phase 11 validity hardening (Tier-1 + WI-7).**
  LOCKED (creating/editing): `vta/eval/baselines.py`, `vta/eval/significance.py`,
  `vta/eval/metrics.py` (append-only), `vta/nodes/rescore.py` (WI-7 contract),
  `vta/nodes/rank.py` (WI-6 GATED weight change: ΔG-primary, LE demoted — supported by
  paired-bootstrap Δ-CI on powered Mpro),
  `scripts/phase11_*.py`, `tests/test_phase11_*.py`, `outputs/phase11/**`.
  Scope this pass: WI-1 (trivial baselines), WI-3 (paired bootstrap + gate redef, remove
  CI-overlap logic — incl. editing Phase-10 `scripts/run_phase10_ensemble.py` +
  `scripts/finalize_phase10.py` to delegate significance to the paired test, and
  regenerating `outputs/phase10/ensemble_benchmark.json` delta from committed ensemble_rows,
  no re-dock), WI-6 (LE-demotion gate), WI-7 (GNINA optional rescore no-op).
  DEFERRED (compute-heavy / tool-absent, not claimed): WI-2 redock RMSD, WI-4 decoy ratio
  ≥30:1 (needs long docking), WI-5 HCV-NI re-benchmark. Coordinate before taking those.

---

## ACTIVE WORK (continued)

- **session-a (fable)** — 2026-06-30 — **Phase 11 WI-4 Job A + WI-5 Job B (extended docks).**
  LOCKED: `scripts/phase11_mpro_30to1.py`, `scripts/phase11_hcv_unmatched.py`,
  `tests/test_phase11_extended.py`, `outputs/phase11/mpro_30to1_*`,
  `outputs/phase11/hcv_ni_unmatched_*`, `vta/data/decoys_cache/hcv_ns5b_ni_unmatched.smi`.
  Job A: 30 actives × 900 real Moonshot inactives (30:1) real-Vina re-dock + baselines/paired.
  Job B: property-unmatched/charge-extrema decoys for HCV-NI (DeepCoy absent) real-Vina
  re-benchmark. SEPARATE outputs — frozen 1:1 headline + locked_benchmark untouched.

## DONE / HANDOFF
- **Phase R reasoning architecture — Stages 1–3 DONE (opus, 2026-07-07).** HemaGuide-inspired
  routing layer: **R1** `dossier_node` (structured target dossier: provenance, pocket
  descriptors incl. metal-in-pocket, benchmarkability), **R2** `triage_router_node` (routes
  full_dock|annotate_only|defer|refuse; downgrades ONLY on positive evidence, so TiLV PB1
  stays full_dock & week2 chain intact; docking.py skips non-full_dock targets, labelled),
  **R4** `verification_node`/`build_verdict` (hard gate: redock<2Å + paired-beats-2D +
  applicability → pass/downgrade/defer/refuse; Mpro → DOWNGRADE), **R5** `annotate_rank_node`
  (labelled 2D-similarity annotation for annotate_only targets, nucleotide caveat, no ΔG
  claim). Report shows a structural verdict banner + annotation section. Wired additively:
  species_resolution→dossier→triage→dock … admet→annotate_rank→verification→report. **259
  tests** (+15). Gates: `outputs/phaseR/gate_R{1,2,R4R5}*.md`. No docking-weight/scorer
  change. **DEFERRED:** R3 PlaybookMemory (stubbed "no_precedent" — needs validated-screen
  corpus + leakage guard), R6 ablation (the proof the workflow is the contribution).
- **Deployment D0 — Scientific readiness gate G1 PASS (opus, 2026-07-03).** New
  `vta/report_envelope.py` — non-removable honesty envelope (disclaimer + pinned frozen
  benchmark hash + per-target grade/CI + trivial-2D-baseline paired verdict + pose
  reliability + nucleotide/metal caveats + out-of-domain flag), rendered as the first report
  card and stored on `state["honesty_envelope"]`. Wired D0.1–D0.4 (baselines/redock/paired/
  LE-demotion) into the output. Fixed two STALE-honesty defects: report notes still said
  ranking was LE-based / conservation a "10% term" (both false post-WI-6) → corrected to
  ΔG-primary, weight-0 annotations; rank.py provenance label fixed. Gate artifact:
  `outputs/deploy/gate_G1_scientific_readiness.md`. **240 tests green.** No ranking-weight
  change; no service (D1+) work started — held for user review after G1 per plan.
- Phase 11 WI-5 reframe (fable) — removed 'IMPOSSIBLE'; HCV-NI finding tied to the purchasable-library scaffold-distinct recipe (property-unmatched/generative/real-nucleotide untested). WI-4 ratio analysis: 30:1 reachable with 30 actives x 900 real inactives (extended dock, not run); DeepCoy absent. 231 tests green.
- Phase 11 WI-2 redocking (fable) — Mpro 7L11/XF1 redocks to 1.654 Å (pose-reliable); triphosphate GTP/CTP RMSD uncomputed (RDKit phosphate valence; spyrmsd absent) — labelled tooling limit; 6Y2E/7K3T apo N/A. 231 tests green.
- Phase 11 Tier-1 (fable session) — WI-1 baselines, WI-3 paired bootstrap + CI-overlap
  purge (incl. Phase-10 scripts), WI-6 LE-demotion (APPLIED: ΔG primary), WI-7 GNINA
  no-op. 228 tests green. DEFERRED: WI-2/4/5 (compute-heavy / DeepCoy absent).


- Phase 9C — Mpro powered-benchmark scaffolding + nucleotide meta-finding (`bce4182`).
- Phase 9D/9E — parent-vs-active-form (NS5B 2XI3 + Mg²⁺) + frozen `locked_benchmark.json` +
  OpenBabel receptor-prep fallback (`54d7505`).
- Phase 10 — ensemble docking on Mpro (7L11 + 6Y2E apo + 7K3T holo); ensemble vs single-
  structure change NOT significant (CIs overlap) (`828cd95`…`97fc0e3`).
- Docs — MASTER_DOCUMENT §15 Validation Results; handoff `docs/claude_code_handoff_2026-06-30.md`
  supersedes 06-29 (`3576b19`).

## NOTES

- Do not re-run finished phases; check `git log --all --oneline` before starting.
- Constraints: no fabrication; no ranking change outside the gate; every metric median+95% CI;
  never pool binding modes.
