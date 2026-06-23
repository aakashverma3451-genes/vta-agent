# VTA-Agent

Autonomous structure-based virtual-screening pipeline (LangGraph): classify → route →
structure → pockets → dock → rescore → rank → admet → report, with an opt-in MD
validation phase (md_select → md_simulate → md_analyze → md_rerank).

## Multi-session collaboration (IMPORTANT)
More than one Claude Code session may run against this repo at once. Sessions are
independent processes and only coordinate through **`.claude/COLLAB.md`** and git.

- **At the start of every task:** read `.claude/COLLAB.md`. Respect any `LOCKED` file
  another session has claimed.
- **Before editing a file:** add an ACTIVE WORK entry claiming it (session id, time,
  file, `LOCKED`).
- **After finishing:** move your entry to DONE / HANDOFF and commit on your branch.
- **Work on your own branch** (`session-a`, `session-b`, …), commit small and often;
  use `git fetch && git log --all --oneline` to see the other session's real changes.

## Commands
- `taxonagent/venv/bin/python -m pytest -q` — run the test suite (use the project venv
  at `../taxonagent/venv`; deps: requests, langgraph, pytest).

## Data sources
External databases the pipeline draws on are catalogued in `vta/data/databases.py`
(machine-readable, with honest INTEGRATED/PLANNED/CATALOGUED status) and mirrored in
`docs/databases.md`. Tool binaries self-discover via `vta/toolconfig.py`. Keep the
registry's `wired_in` flags honest: only mark INTEGRATED what a node actually reads.

## Conventions
- Every external call goes through an injectable network seam (so tests stay offline)
  with a committed cache or labelled fallback — never crash, degrade honestly.
- Tests are hermetic (monkeypatch the seams); add a test with each new node/branch.
