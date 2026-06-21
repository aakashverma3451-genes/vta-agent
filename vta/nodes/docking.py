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

from vta.state import VTAState

# Stand-in ligand set. The 5 named compounds are real nucleoside-analog RdRp
# inhibitors used as POSITIVE CONTROLS — in Phase 2 (real docking) they must rank
# in the top decile, or the pipeline is mis-calibrated. The rest are placeholders
# for a real ChEMBL antiviral query.
POSITIVE_CONTROLS = ["remdesivir", "favipiravir", "ribavirin", "sofosbuvir", "molnupiravir"]
MOCK_LIGANDS = POSITIVE_CONTROLS + [f"chembl_{i:04d}" for i in range(95)]   # 100 total

_TOP_POCKETS = 3


def _seed(*parts: str) -> int:
    h = hashlib.md5("|".join(parts).encode()).hexdigest()
    return int(h[:8], 16)


def docking_node_mock(state: VTAState) -> VTAState:
    results = []
    for protein, plist in (state.get("pockets") or {}).items():
        for pocket in plist[:_TOP_POCKETS]:
            for lig in MOCK_LIGANDS:
                rng = random.Random(_seed(state["run_id"], protein, str(pocket["id"]), lig))
                dG = round(rng.uniform(-9.0, -5.0), 2)
                rmsd = round(rng.uniform(0.5, 3.0), 2)
                heavy_atoms = rng.randint(20, 45)
                le = round(dG / heavy_atoms, 3)
                results.append({
                    "protein": protein,
                    "pocket": pocket["id"],
                    "ligand": lig,
                    "dG": dG,
                    "rmsd": rmsd,
                    "le": le,
                    "heavy_atoms": heavy_atoms,
                    "conservation": pocket["conservation"],
                })
    state["docking_results"] = results
    state["audit_trail"].append(f"[MOCK] Docking: {len(results)} ligand-pocket pairs")
    state["versions"]["docking"] = "MOCK-vina-shaped"
    return state
