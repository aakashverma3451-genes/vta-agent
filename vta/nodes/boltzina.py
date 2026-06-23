"""boltzina_node — Boltzina DL affinity prediction seam (§1.2 / arXiv:2508.17555).

Boltzina takes Vina-docked poses and re-scores them with a lightweight Boltz-2
variant, yielding a binding-likelihood score that outperformed both Vina and GNINA
on the MF-PCBA benchmark. Where GNINA rescores PDBQT poses via a binary, Boltzina
is a Python package that takes protein PDB + ligand SMILES — so it can annotate
even mock-docking records (which carry real SMILES but no saved pose files).

Input per docking record:  protein PDB path (from state["structures"])  +  SMILES
Output:                    boltzina_score  (0–1, higher = more likely to bind)

ANNOTATION-ONLY by design (same discipline as rescore_node / admet_node).
The term is not yet wired into rank.py — blending requires a gate recalibration
on a GPU host. Surface as a column/flag first, weight deliberately later.

Auto-detects the `boltzina` package; degrades to labelled skip when absent.

NOTE ON API: Boltzina (arXiv:2508.17555, Aug 2025) was published after the
training cutoff. The `_boltzina_predict()` stub below reflects the most likely
API from the paper description; verify against the released package and update
the call if the signature differs.
"""
from __future__ import annotations

from vta.state import VTAState


def _boltzina_available() -> bool:
    try:
        import boltzina  # noqa: F401
        return True
    except Exception:
        return False


def _boltzina_predict(protein_pdb: str, smiles: str) -> float | None:
    """Call Boltzina to predict binding likelihood (0–1) for one protein-ligand pair.

    Verify the exact API against the released package; update if needed.
    The paper describes: protein structure + ligand SMILES → binding probability.
    """
    try:
        from boltzina import BoltzinaModel
        model = BoltzinaModel()
        result = model.predict(protein_pdb=protein_pdb, ligand_smiles=smiles)
        # result may be a float, a dict, or a dataframe depending on the release.
        if isinstance(result, (int, float)):
            return float(result)
        if isinstance(result, dict):
            return float(result.get("score") or result.get("binding_probability", 0))
        return None
    except Exception:
        return None


def boltzina_node(state: VTAState) -> VTAState:
    """Annotate docking records with Boltzina binding-likelihood scores."""
    rows = state.get("docking_results") or []
    if not rows:
        return state

    if not _boltzina_available():
        state["audit_trail"].append(
            "[skip] Boltzina: package not installed (pip install boltzina)")
        state["versions"]["boltzina"] = "skipped"
        return state

    structures = state.get("structures") or {}
    scored, skipped = 0, 0

    for r in rows:
        protein = r.get("protein")
        smiles = r.get("smiles")
        pdb_path = (structures.get(protein) or {}).get("pdb_path")

        if not (smiles and pdb_path):
            skipped += 1
            continue

        score = _boltzina_predict(pdb_path, smiles)
        if score is not None:
            r["boltzina_score"] = round(score, 4)
            scored += 1
        else:
            skipped += 1

    if scored:
        best = max((r for r in rows if "boltzina_score" in r),
                   key=lambda r: r["boltzina_score"])
        state["versions"]["boltzina"] = "Boltzina (Boltz-2 variant, arXiv:2508.17555)"
        state["audit_trail"].append(
            f"Boltzina: scored {scored} records (annotation-only, ranking unchanged); "
            f"best score {best['boltzina_score']} ({best.get('ligand')})"
        )
    else:
        state["versions"]["boltzina"] = "skipped (no scoreable records)"
        state["audit_trail"].append(
            "[skip] Boltzina: package present but no records had protein PDB + SMILES")

    return state
