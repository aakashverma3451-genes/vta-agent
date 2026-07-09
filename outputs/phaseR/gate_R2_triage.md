# Phase R — Stage 2 gate (R2 TriageRouter): PASS

**Date:** 2026-07-07 · **Owner:** session-a (opus)

**Go criterion (from plan):** the router correctly downgrades the known-hard cases in a
labelled pilot; the central lesson (route the target, don't dock everything) is wired.

- `vta/nodes/triage.py` (`triage_router_node`) reads the R1 dossier + R3 playbook stub and
  writes `triage_decision` per protein ∈ {full_dock | annotate_only | defer | refuse} with a
  rationale. Wired additively as `dossier → triage → dock`; `docking.py` skips (labelled) any
  protein the router left off `full_dock`. R3 PlaybookMemory is stubbed to "no_precedent".
- Verified routing (`tests/test_triage.py`):
  - **Mpro** → `full_dock` (experimental, in-domain, powered).
  - **HCV NS5B NI** → `annotate_only` via `benchmarkability=un_benchmarkable` (the reliable
    trigger — matched-decoy validation is infeasible for this class).
  - **TiLV PB1** → `full_dock` (metal-dependent polymerase biology, but only *underpowered*
    and no metal resolved in its pocket) → the week2 e2e chain keeps its 20 docked leads.
  - **Low binding-site pLDDT predicted** → `defer` (60) / `refuse` (40).
  - **Metal in the detected pocket** or an explicit **nucleotide/charged/covalent screened
    ligand class** → `annotate_only`.
- Invariant held: the router changes WHICH method runs, never the docking weights/scoring.
  Default is `full_dock`; downgrades fire only on positive, run-available evidence.

**Evidence:** 251 tests pass (+7). R4 verification gate + R5 annotation ranker next (Stage 3);
report at end of Stage 3 before the R3 playbook corpus.
