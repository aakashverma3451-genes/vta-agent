# Scientific Validation Status

## Completed

| Validation Item | Status | Evidence |
|-----------------|--------|----------|
| Test suite | Completed | `143 passed` in this workspace before environment refactor |
| Positive control benchmark | Completed | User-reported PASS, 3/4 controls in Top 5 |
| AutoDock Vina integration | Completed | User-reported AutoDock Vina 1.2.5 validated |
| FPocket integration | Completed | User-reported FPocket validated |
| FEP seam tests | Completed | Hermetic tests for skipped and stubbed FEP execution |

## Pending Real-Tool / Production Validation

| Validation Item | Status | Blocker |
|-----------------|--------|---------|
| ADMET-AI | Pending | Install `admet-ai` and run real report |
| GNINA | Pending | GPU host and `gnina` binary |
| Boltzina | Pending | Package/API verification and GPU-capable environment |
| Molecular dynamics | Pending | OpenMM/MDAnalysis environment and compute time |
| FEP/ABFE | Pending | `openfe-abfe` runner and long compute |

## Rule For Status Changes

Move an item to completed only when:

- The exact command is recorded.
- Tool versions are recorded.
- Output artifacts are saved or reproducible.
- The report or validation JSON reflects the result.

