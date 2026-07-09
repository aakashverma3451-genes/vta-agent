# Phase R — Stage 1 gate (R1 DossierBuilder): PASS

**Date:** 2026-07-07 · **Owner:** session-a (opus)

**Go criterion (from plan):** the dossier is correct on the three wired targets.

- `vta/nodes/dossier.py` (`dossier_node`) assembles one `target_dossier` per protein from the
  structure/pocket nodes + a curated target-knowledge table, inserted at
  `species_resolution → dossier → dock` (additive; no edges removed).
- Verified on the three wired targets (`tests/test_dossier.py`):
  - **Mpro** → protease, `metal_dependence=False`, benchmarkability **powered**.
  - **HCV NS5B** → polymerase, `metal_dependence=True`, benchmarkability **un_benchmarkable**
    (from the frozen benchmark grade — this is the reliable downgrade trigger for R2).
  - **TiLV PB1** → metal-dependent polymerase biology but benchmarkability **underpowered**
    (NOT un_benchmarkable) → stays dockable, so the week2 e2e chain keeps its 20 docked leads.
- Honesty: unknown fields are `null` with a reason, never fabricated; `binding_site_plddt`
  carries the mean-pLDDT proxy with an explicit note (per-binding-site pLDDT not wired;
  Scardino 2023 — global pLDDT necessary-not-sufficient).

**Evidence:** 244 tests pass (+4). No docking-weight/scorer change. R2 next.
