# COLLAB — shared whiteboard for parallel Claude Code sessions

Two (or more) `claude` sessions are independent processes; they only see each other
through this file and through git. **This file is the live status board.**

## Protocol (every session, every task)
1. **Before** starting a task: read this file. If another session has LOCKED a file
   you need, pick something else or coordinate.
2. **Claim** your work: add a line under ACTIVE WORK with your session id, time,
   the file(s) you'll touch, and `LOCKED` for any file you're editing.
3. **After** finishing: move the entry to DONE / HANDOFF and `git commit` on your
   branch so the change is real, not just announced here.
4. Keep ACTIVE WORK to what's truly in flight. Stale locks block your partner.

## Branch discipline — git worktrees (each session = own folder)
Two terminals in the SAME folder share one working tree and one HEAD, so plain
`git checkout` can't isolate them. Each session instead has its own **worktree**:

- **Session A** → `/Users/dna/Desktop/vtaagent/vta-agent`            branch `session-a`
- **Session B** → `/Users/dna/Desktop/vtaagent/vta-agent-session-b`  branch `session-b`

Session B's terminal must `cd` into the `-session-b` folder. Commit small and often;
merge to `main` when a unit is done. Use `git log --all --oneline -15` to see the
other session's real commits — git is the source of truth; this file is just intent.

A non-blocking pre-commit hook reminds you to update this file (it never blocks).

---

## ACTIVE WORK
- _(none — claim your task here)_

## SPECS — FOR ANY IMPLEMENTER (posted by the research lead session)
> Researched, buildable. CLAIM in ACTIVE WORK before starting. Ping the research-lead
> session via this file if the validation gate moves — weight re-tuning is mine to call.

### SPEC #1 — real conservation scoring (retire the `0.5` placeholder)  [UNCLAIMED]
**Goal.** Replace `pockets.py`'s `_CONSERVATION_PLACEHOLDER = 0.5` with a real per-pocket
score in [0,1] so `rank.py`'s EXISTING 10% conservation term carries signal. Higher =
pocket residues are evolutionarily conserved across viral homologs (harder to escape by
mutation — a genuine antiviral druggability axis). This is the most visible fake in the
gate output (`cons 0.5` for every ligand).

**Metric (grounded).** Jensen–Shannon divergence per MSA column vs a background AA
distribution — Capra & Singh 2007, *Bioinformatics* 23(15):1875, doi:10.1093/
bioinformatics/btm270. Best-performing conservation measure for functional residues AND
natively bounded [0,1] (no renorm). Ref impl: github.com/shin-kinos/cons-capra07.

**Design.** New node `vta/nodes/conservation.py`, wired `pockets → conservation → dock`
(docking.py ALREADY copies `pocket["conservation"]`; do NOT touch docking/rank logic).
Per protein with a structure:
  1. **MSA seam (injectable):** `fetch_homolog_msa(sequence, taxon) -> list[str]` →
     gapped homolog rows, row 0 = target. Real path: homolog search (NCBI/UniProt) +
     MAFFT; COMMIT a TiLV-PB1 MSA cache so the gate is reproducible offline. No MSA →
     keep 0.5 + audit `[skip] conservation: no MSA` (graceful-fallback discipline).
  2. **JSD (pure, dependency-free):** `jensen_shannon_conservation(cols) -> list[float]`
     in [0,1]; hardcode standard AA background freqs. (Capra neighbour-window = v2.)
  3. **Pocket residues by PROXIMITY (universal):** the experimental site is only a
     `center`, no residue list — so take residues with any atom within R=6 Å of
     `pocket["center"]` in the structure PDB. Map structure resid → target-seq index →
     MSA column. ⚠ NUMBERING HAZARD: experimental chains (8PSO) aren't guaranteed
     1-indexed/contiguous — align observed residues to the provided sequence, don't
     assume an offset. Add a test for this.
  4. **Aggregate:** pocket conservation = mean JSD over pocket residues (round [0,1]);
     overwrite `pocket["conservation"]`.

**Files.** +`vta/nodes/conservation.py`, +`tests/test_conservation.py`, edit
`vta/graph.py` (insert node, both phases), +committed MSA cache under `vta/data/`,
`conftest.py` monkeypatch `fetch_homolog_msa -> None` (keep CI hermetic), register the
MSA source in `vta/data/databases.py`, tweak `pockets.py` + `report.py` placeholder note.

**Tests (hermetic).** JSD on toy MSA (identical col → ~1.0, uniform → ~0); proximity
mapping on toy PDB incl. the numbering-offset case; seam-absent → 0.5 + label; node
end-to-end with stubbed MSA → non-0.5 on the pocket.

**ACCEPTANCE (MANDATORY — this is a ranking term).** Full suite green + new tests; then
RE-RUN `scripts/validate_controls.py` — gate MUST still PASS (≥3/4 controls top-5) and
paste the new control table into your DONE entry. If real conservation breaks the gate,
DO NOT ship the regression — flag the research lead; `rank.py WEIGHTS["conservation"]
= 0.10` is the recalibration knob. The report's "all 0.5 → planned" note must disappear
on the real path.

### SPEC #2 — prmtop handoff: make MM-GBSA runnable (md_simulate emits an Amber topology)  [UNCLAIMED]
**Goal.** `md_analyze`'s MM-GBSA is wired and reads `result["parm"]`, but `md_simulate`
never emits one — so production MM-GBSA can't run. Close the handoff: have `md_simulate`
write a complex Amber topology and add its path to each completed result.

**Tool (grounded).** ParmEd bridges OpenMM→Amber: `parmed.openmm.load_topology(topology,
system, xyz=positions)` → `Structure`, then `.save("complex.prmtop")` + `.save(
"complex.inpcrd")`. ParmEd natively translates an OpenMM System to Amber prmtop/inpcrd
(parmed.github.io/ParmEd → amber package). Needs the PARAMETERISED system (forces set).

**Design — `vta/nodes/md_simulate.py` ONLY.** Add `_check_parmed()` seam. After the
production run, for each completed sim: ParmEd-save `<run_dir>/complex.prmtop` (+inpcrd),
set `result["parm"] = <prmtop path>`. Graceful: ParmEd absent OR save throws → omit
`parm` (md_analyze already prints the honest "no prmtop" note — DON'T duplicate it).
Keep the existing OpenMM-absent skip untouched. Do NOT edit md_analyze.

**Files.** `vta/nodes/md_simulate.py`, tests in `tests/test_md.py`. **Disjoint from
SPEC #1 and #3** — no shared files. (md_simulate isn't in the default graph, so no
conftest/graph edits.)

**Tests (hermetic).** Monkeypatch the system build + `_check_parmed`/ParmEd save:
assert a completed `result` carries `parm` when ParmEd present; omitted when absent.
Then feed that result to `md_analyze._compute_mmgbsa` with a stubbed `_run_mmgbsa_tool`
→ confirm it no longer returns the "no Amber topology" note.

**Acceptance.** Full suite green + new tests. No validation-gate impact (MD is opt-in,
not in the default-graph gate). Note in DONE that production MM-GBSA on a GPU host now
has its topology.

### SPEC #3 — refresh the ARCHITECTURE diagram to match the real graph  [UNCLAIMED]
**Goal.** `ARCHITECTURE.md` A.2 state-machine (and A.1/A.3) are stale: they show
`structure → pockets → dock → rank`, label `pockets`/`dock` as MOCK (both are REAL
auto-discovered tools now), and omit `proteinttt`, `rescore`, `boltzina`, `admet`, the
MD phase, and (pending #1) `conservation`. A half-fix is worse than honest-stale, so do
the WHOLE pass.

**Design — `ARCHITECTURE.md` ONLY (docs, no code/tests).** Make the diagram nodes+edges
EXACTLY match `build_app()` in `vta/graph.py` HEAD: classify → router → structure →
proteinttt → pockets → dock → rescore → boltzina → rank → admet → report, plus the
opt-in MD phase (md_select → md_simulate → md_analyze → md_rerank). Update the
green/orange legend: drop "mock" for pockets/dock; represent the honest real-vs-labelled-
fallback reality instead. Refresh A.3 state-growth to include the fields added since
(structures/proteinttt, rescore/boltzina annotations, admet, md_*).

**Files.** `ARCHITECTURE.md` only. **Disjoint from #1 and #2.**

**Sequencing.** Soft-depends on SPEC #1: if conservation has merged, include the
`conservation` node; if not, add it with a "(pending SPEC #1)" note and finalise once #1
lands. Acceptance: every node/edge in the diagram exists in `graph.py`; no node labelled
mock that auto-discovers a real tool.

## DONE / HANDOFF
- **[DB task FINISHED — session-a]** New `vta/data/pubchem.py` (PUG-REST name→SMILES);
  wired as a fallback in `ligands.py` so a ChEMBL miss is recovered from PubChem, not
  dropped. SwissADME reclassified CATALOGUED (web form, no public API; ADMET is done
  locally via admet_ai). Registry now 6 INTEGRATED / 0 PLANNED / 23 CATALOGUED — nothing
  dangling. +6 pubchem tests, +2 registry tests; full suite 82 passed (7 pre-existing
  classify_genome failures unrelated).
- **[MM-GBSA — md_analyze]** Full MM-GBSA binding free energy (ΔTOTAL) in
  `vta/nodes/md_analyze.py`: gmx_MMPBSA→MMPBSA.py auto-detect, real
  FINAL_RESULTS_MMPBSA.dat parser, subprocess seam, graceful labelled skip when the
  engine OR the Amber topology (prmtop) is absent. Annotation-only (NOT in md_rerank).
  HANDOFF: md_simulate should emit a `parm`/`prmtop` (ParmEd) so production runs have
  a topology — node already reads `result["parm"]`. +8 tests in `tests/test_md.py`
  (parser + skip/missing-topology/compute/failure + node-level). Full suite 82 passed.
- **[UniProt resolver — session-a]** New `vta/data/uniprot.py` (accession resolution +
  sequence fetch, seam+cache+offline fallback). Opt-in hook in `structure.py`: a
  `uniprot_query` on a >400aa orphan protein resolves an accession → AlphaFold DB.
  Registry marks `uniprot` INTEGRATED. +7 uniprot tests, +2 structure tests; full
  suite 67 passed (7 pre-existing classify_genome failures unrelated). Committing now.
- **[setup]** Created this whiteboard + project `CLAUDE.md` check-in rule.
- **[structure stage]** AlphaFold DB fallback wired into `vta/nodes/structure.py`
  (>400aa proteins with a UniProt accession now fetch from AlphaFold DB ahead of
  local Boltz-2; graceful 404 degrade). Registry `vta/data/databases.py` marks
  `alphafold` INTEGRATED. Tests: `tests/test_structure.py` (+2), all 19 pass.
- **[databases]** Added `docs/databases.md`, `vta/data/databases.py` registry, and
  `tests/test_databases.py` (7 pass). Committed in `cac9492`.
- **[structure merge — session-a]** RESOLVED the contended `structure.py`: kept the
  full cascade ESMFold (≤400aa) → AlphaFold DB (>400aa w/ accession) → Boltz-2
  (>400aa local) → honest refuse. Added 4 Boltz-2 tests; `tests/test_structure.py`
  now 12/12, full suite 59/59. Committed on `session-a`.
- **[§2.3 ProteinTTT — session-a]** New `vta/nodes/proteinttt.py`: refines
  low-pLDDT (<70) ESMFold folds via test-time training; experimental/AlphaFold/
  Boltz-2 left untouched; graceful skip w/o the `proteinttt` package (this box).
  Wired `structure → proteinttt → pockets`. +6 tests; full suite 65/65. Committed
  on `session-a`.

## CONTENDED FILES (heads up)
- _(none — `structure.py` merge resolved; both cascades committed)_

## BLOCKED / NEEDS DECISION
- _(none)_
