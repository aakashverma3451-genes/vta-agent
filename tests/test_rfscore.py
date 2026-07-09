"""Hermetic tests for the in-house RF-Score rescorer (Phase S). Tiny synthetic structures."""
from __future__ import annotations

import numpy as np

from vta.eval.rfscore import (FEATURE_NAMES, murcko_scaffold, read_atoms,
                              rfscore_features, scaffold_cv_predict)

_RECEPTOR_PDB = (
    "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
    "ATOM      2  N   ALA A   1       0.000   0.000  20.000  1.00  0.00           N\n"
    "END\n")
# ligand PDBQT: two carbons near the receptor CA + a 'G0' macrocycle pseudo-atom (must be skipped)
_LIGAND_PDBQT = (
    "ATOM      1  C   UNL     1       0.000   0.000   0.000  1.00  0.00     0.000 C\n"
    "ATOM      2  C   UNL     1       1.000   0.000   0.000  1.00  0.00     0.000 C\n"
    "ATOM      3  G0  UNL     1       5.000   0.000   0.000  1.00  0.00     0.000 G0\n"
    "ENDMDL\n")


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def test_read_atoms_skips_pseudo_and_hydrogens(tmp_path):
    lig = read_atoms(_write(tmp_path, "l.pdbqt", _LIGAND_PDBQT), "pdbqt")
    assert set(lig) == {"C"} and lig["C"].shape == (2, 3)   # 2 carbons, G0 dropped
    rec = read_atoms(_write(tmp_path, "r.pdb", _RECEPTOR_PDB), "pdb")
    assert rec["C"].shape == (1, 3) and rec["N"].shape == (1, 3)


def test_rfscore_features_counts_contacts_within_cutoff(tmp_path):
    rec = read_atoms(_write(tmp_path, "r.pdb", _RECEPTOR_PDB), "pdb")
    f = rfscore_features(rec, _write(tmp_path, "l.pdbqt", _LIGAND_PDBQT))
    feat = dict(zip(FEATURE_NAMES, f))
    # protein C at origin ↔ both ligand C within 12 Å → 2 contacts; protein N is 20 Å away → 0
    assert feat["C-C"] == 2.0
    assert feat["N-C"] == 0.0
    assert f.shape == (36,)


def test_rfscore_features_none_on_unreadable(tmp_path):
    rec = read_atoms(_write(tmp_path, "r.pdb", _RECEPTOR_PDB), "pdb")
    assert rfscore_features(rec, str(tmp_path / "missing.pdbqt")) is None


def test_murcko_scaffold_groups_analogs_together():
    # two biphenyl analogs (differ only in an alkyl tail) share a generic scaffold; benzene doesn't
    s1 = murcko_scaffold("c1ccc(-c2ccccc2)cc1CC")
    s2 = murcko_scaffold("c1ccc(-c2ccccc2)cc1CCC")
    assert s1 and s1 == s2
    assert murcko_scaffold("c1ccccc1") != s1        # single-ring scaffold ≠ biphenyl scaffold


def test_scaffold_cv_predict_deterministic_and_covers_all():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(30, 36))
    y = X[:, 0] * 2 + rng.normal(scale=0.1, size=30)   # learnable signal
    groups = [f"g{i % 6}" for i in range(30)]          # 6 scaffold groups
    a = scaffold_cv_predict(X, y, groups, n_splits=3, seed=0)
    b = scaffold_cv_predict(X, y, groups, n_splits=3, seed=0)
    assert np.allclose(a, b)                            # deterministic under fixed seed
    assert not np.isnan(a).any()                        # every row predicted out-of-fold
    # a signal-following model should correlate with truth out-of-fold
    assert np.corrcoef(a, y)[0, 1] > 0.3
