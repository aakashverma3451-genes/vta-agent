# Phase R — Stage 3 gate (R4 VerificationNode + R5 annotation ranker): PASS

**Date:** 2026-07-07 · **Owner:** session-a (opus)
**This is the reporting checkpoint** — per the plan, report at the end of Stage 3 before
building the R3 playbook corpus.

**Go criteria (from plan):**
1. Mpro emits **DOWNGRADE** (not a docking ranking).
2. **Zero** docking-enriched claims escape a failed paired-baseline gate.
3. `annotate_only` targets return a labelled ligand-based annotation, never a docking-precise
   enrichment number.

## What was built
- `vta/nodes/verification.py` (`verification_node` / pure `build_verdict`) — the hard gate on
  the single edge into reporting. Three deterministic gates: **redock-RMSD < 2 Å**
  (redock_validation.json), **paired-baseline-beats-2D** (baselines_and_gate.json
  `structure_based_enrichment_demonstrated`, the D0.3 gate promoted to a hard pre-report gate),
  and **applicability-domain** (must be `full_dock`). Run verdict = worst-case across
  docking-claim targets: `pass | downgrade | defer | refuse | not_applicable`.
- `vta/nodes/annotate_rank.py` (`annotate_rank_node`) — the `annotate_only` executor: ranks a
  supplied known-actives library by 2-D ECFP4 similarity, explicitly tagged
  `ligand_based_2d_similarity_annotation`; else a labelled `annotation_only` record. Carries
  the Phase-9D nucleotide/metal caveat. Never attaches a ΔG/enrichment number.
- `vta/report.py` — structural verdict **banner** (downgrade in red) + a **ligand-based
  annotations** section. Wired additively: `admet → annotate_rank → verification → report`
  (and `md_rerank`/`fep → verification → report`).

## Evidence the criteria are met
- **(1) & (2)** `build_verdict` on an Mpro docking-lead state → `verdict = downgrade`,
  `beats_2d_baseline = False`; an *untested* target (no baseline artifact) also downgrades
  (conservative — an unverifiable enrichment claim is never passed). The report renders
  "⚠ structure-based enrichment NOT demonstrated — docking ranking DOWNGRADED"
  (`tests/test_verification.py`, verified against the real committed artifacts).
- **(3)** End-to-end offline drive: HCV NS5B (un_benchmarkable) → triage `annotate_only` →
  **docking skipped** (0 results, labelled audit line) → `ligand_based_annotation` with the
  nucleotide caveat → report shows the annotation section, not a docking ranking.
- Non-regression: **259 tests pass** (+8 this stage). No docking-weight/scorer change; the
  router/verifier change only WHICH method runs and whether a claim may be made.

## Honest limitations carried forward
- R5's 2D-similarity annotation needs a supplied known-actives library; without one it returns
  a labelled "no validated ranking" record (honest, but not a ranking).
- The paired-baseline gate is per-target and only as current as `baselines_and_gate.json`
  (a stale benchmark would need the D10 revalidation path, not yet built).
- An LLM chain-of-verification layer over the narrative is designed but not implemented; by
  design it may never override these physics/statistics gates.

## Deferred (not in this session, per plan)
- **R3 PlaybookMemory** — stubbed "no_precedent"; needs a validated-screen corpus + the
  leakage guard before it pays off.
- **R6 component ablation** — the L0→full-agent study that would *prove* the workflow (not the
  docking) is the contribution. This remains a hypothesis until R6 runs.
