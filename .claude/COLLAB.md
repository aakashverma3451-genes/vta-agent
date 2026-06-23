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

## DONE / HANDOFF
- **[setup]** Created this whiteboard + project `CLAUDE.md` check-in rule.
- **[structure stage]** AlphaFold DB fallback wired into `vta/nodes/structure.py`
  (>400aa proteins with a UniProt accession now fetch from AlphaFold DB ahead of
  local Boltz-2; graceful 404 degrade). Registry `vta/data/databases.py` marks
  `alphafold` INTEGRATED. Tests: `tests/test_structure.py` (+2), all 19 pass.
- **[databases]** Added `docs/databases.md`, `vta/data/databases.py` registry, and
  `tests/test_databases.py` (7 pass). Uncommitted on `master` as of this writing.

## CONTENDED FILES (heads up)
- `vta/nodes/structure.py` — Session A added the AlphaFold-DB fallback; a Boltz-2
  branch from the other terminal landed in the same file. Both are interleaved and
  uncommitted. Whoever commits it should keep BOTH cascades (ESMFold → AlphaFold →
  Boltz-2) and re-run `tests/test_structure.py`.

## BLOCKED / NEEDS DECISION
- _(none)_
