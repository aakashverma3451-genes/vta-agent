# Phase 11 WI-2 — Redocking / pose-reproduction validation

RMSD: RDKit GetBestRMS (symmetry-corrected, heavy-atom); spyrmsd not installed. Flag threshold: 2.0 A.

| Target | PDB | Ligand | Receptor prep | Redock dG | RMSD (A) | Verdict |
|---|---|---|---|---|---|---|
| Mpro_7L11 | 7L11 | XF1 | openbabel_fallback | -9.317 | 1.654 | pose-reliable |
| HCV_NS5B_2XI3 | 2XI3 | GTP | openbabel_fallback | -8.5 | None | RMSD computation failed |
| TiLV_8PSO | None | — | — | — | None | crystal ligand bond-order assignment failed |
| Mpro_6Y2E | 6Y2E | — | — | — | None | N/A (apo) |
| Mpro_7K3T | 7K3T | — | — | — | None | N/A (apo) |

## Interpretation

The powered-benchmark structure Mpro 7L11 redocks its native non-covalent ligand XF1 to 1.654 A (< 2.0 A) — the production Vina/Meeko+OpenBabel protocol reproduces the crystallographic pose, so the WI-1 result (Vina not beating 2D-similarity) is NOT a gross pose-search failure; it is a ranking/benchmark limitation. RMSD for the triphosphate ligands (2XI3/GTP, 8PSO/CTP) was NOT computed: RDKit bond-order assignment fails on the triphosphate (phosphate valence) and spyrmsd is not installed — a TOOLING limitation on those two ligands, not evidence of an unreliable pose. 6Y2E and 7K3T(chain A) are apo.

Any target with RMSD > 2.0 A is flagged pose-unreliable (Phase 11 WI-2).
