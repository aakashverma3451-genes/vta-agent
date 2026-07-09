"""Hermetic tests for WI-2 redocking helpers (no network, no docking)."""
from __future__ import annotations

import copy

from rdkit import Chem
from rdkit.Chem import AllChem

from scripts.phase11_redock import heavy_atom_rmsd, redock_target


def _embedded(smiles: str):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, AllChem.ETKDGv3())
    return Chem.RemoveHs(m)


def test_heavy_atom_rmsd_self_is_zero():
    m = _embedded("CCOc1ccccc1C#N")
    assert heavy_atom_rmsd(m, copy.deepcopy(m)) == 0.0


def test_heavy_atom_rmsd_realigns_rigid_translation():
    m = _embedded("CCOc1ccccc1C#N")
    m2 = copy.deepcopy(m)
    conf = m2.GetConformer()
    for i in range(m2.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (p.x + 5.0, p.y, p.z))
    # GetBestRMS superposes first → rigid translation gives ~0 RMSD.
    assert heavy_atom_rmsd(m, m2) < 0.01


def test_redock_target_apo_is_na_without_docking():
    # Apo branch returns before any network/docking, so vina bin is irrelevant.
    out = redock_target({"name": "Mpro_6Y2E", "pdb": "6Y2E", "apo": True,
                         "note": "apo Mpro"}, vina="/nonexistent/vina")
    assert out["status"] == "N/A (apo)"
    assert out["rmsd"] is None and out["pose_reliable"] is None
