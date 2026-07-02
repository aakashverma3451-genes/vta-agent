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

- **session-a (fable)** — 2026-06-30 — **Phase 11 WI-2 redocking validation.** LOCKED:
  `scripts/phase11_redock.py`, `tests/test_phase11_redock.py`, `outputs/phase11/redock_*`.
  Redocks native co-crystal ligands (7L11/XF1, 2XI3/GTP, 8PSO/CTP) with the production
  Vina protocol; symmetry-corrected heavy-atom RMSD via RDKit GetBestRMS (spyrmsd absent);
  6Y2E + 7K3T flagged apo/N-A. Diagnoses the WI-1 finding (does Vina reproduce known poses?).

## DONE / HANDOFF
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
