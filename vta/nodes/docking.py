"""docking_node (MOCK) — virtual screening, AutoDock-Vina-shaped output.

Plan Task 2.3 (the docking half). Mocked in Phase 1; real Vina + RDKit/Meeko ligand
prep + ChEMBL queries in Phase 2. The fabricated records match real Vina output:
binding free energy `dG` (~ -5 to -9 kcal/mol for plausible binders), pose `rmsd`,
and ligand efficiency `le = dG / heavy_atoms`. The ranking node consumes this exact
shape, so swapping in real Vina changes nothing downstream.

Determinism: scores are seeded from (run_id, protein, pocket, ligand), so a run is
fully reproducible — the same auditability requirement that makes this whole
pipeline defensible.
"""
from __future__ import annotations

import hashlib
import random

from vta.data.ligands import load_ligands
from vta.state import VTAState

_TOP_POCKETS = 3


def _seed(*parts: str) -> int:
    h = hashlib.md5("|".join(parts).encode()).hexdigest()
    return int(h[:8], 16)


def docking_node_mock(state: VTAState) -> VTAState:
    """Mock docking over the REAL ligand library (Phase-2 ChEMBL data), mock scores.

    Ligand IDENTITY is now real (ChEMBL id + name + SMILES from `load_ligands`), so
    Phase 2 only has to replace the fabricated dG/rmsd with real AutoDock Vina output
    — the records already carry the SMILES that real ligand prep needs.
    """
    ligands = load_ligands()
    results = []
    for protein, plist in (state.get("pockets") or {}).items():
        for pocket in plist[:_TOP_POCKETS]:
            for lig in ligands:
                rng = random.Random(
                    _seed(state["run_id"], protein, str(pocket["id"]), lig["chembl_id"]))
                dG = round(rng.uniform(-9.0, -5.0), 2)
                rmsd = round(rng.uniform(0.5, 3.0), 2)
                heavy_atoms = rng.randint(20, 45)
                le = round(dG / heavy_atoms, 3)
                results.append({
                    "protein": protein,
                    "pocket": pocket["id"],
                    "ligand": lig["name"],
                    "ligand_id": lig["chembl_id"],
                    "smiles": lig["smiles"],
                    "positive_control": lig["positive_control"],
                    "dG": dG,
                    "rmsd": rmsd,
                    "le": le,
                    "heavy_atoms": heavy_atoms,
                    "conservation": pocket["conservation"],
                })
    state["docking_results"] = results
    state["audit_trail"].append(
        f"[MOCK] Docking: {len(results)} pairs over {len(ligands)} real ligands")
    state["versions"]["docking"] = "MOCK-vina-shaped over ChEMBL ligands"
    state["versions"]["ligands"] = f"ChEMBL x{len(ligands)}"
    return state
