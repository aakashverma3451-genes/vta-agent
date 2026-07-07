"""Hermetic tests for the activity-cliff scoring (Phase S). Synthetic compounds; RDKit only."""
from __future__ import annotations

import math

from vta.eval.cliffs import (
    bootstrap_paired,
    find_cliff_pairs,
    pic50_from_um,
    score_cliff_pairs,
)


def test_pic50_from_um():
    assert pic50_from_um(1.0) == 6.0                 # 1 µM → pIC50 6
    assert abs(pic50_from_um(0.1) - 7.0) < 1e-9      # 100 nM → pIC50 7
    assert pic50_from_um(None) is None
    assert pic50_from_um(0) is None                  # non-positive → labelled skip


def _cpds():
    # two identical-SMILES compounds (Tanimoto 1.0) with a large potency gap = a cliff pair;
    # a third, dissimilar compound provides the kNN neighbour.
    return [
        {"id": "A_potent", "smiles": "CCO", "pic50": 8.0, "dG": -9.0},
        {"id": "A_weak",   "smiles": "CCO", "pic50": 5.0, "dG": -6.0},
        {"id": "other",    "smiles": "c1ccccc1", "pic50": 6.0, "dG": -7.0},
    ]


def test_find_cliff_pairs_flags_similar_high_delta_only():
    pairs, fps = find_cliff_pairs(_cpds(), sim_threshold=0.7, dpic50_threshold=1.0)
    assert len(pairs) == 1
    i, j, sim, dp = pairs[0]
    assert {i, j} == {0, 1} and sim == 1.0 and dp == 3.0    # the identical-SMILES cliff


def test_score_cliff_pairs_vina_correct_when_more_potent_binds_stronger():
    cpds = _cpds()
    pairs, fps = find_cliff_pairs(cpds, sim_threshold=0.7, dpic50_threshold=1.0)
    scored = score_cliff_pairs(pairs, cpds, fps, k=3)
    # A_potent (pIC50 8) has the lower ΔG (-9) → Vina ranks the cliff correctly
    assert scored["vina"] == [1]
    # 2D-kNN (excluding both pair members) can only see the dissimilar 'other' → cannot resolve
    assert scored["twod"] == [0]


def test_score_drops_pairs_missing_dG():
    cpds = _cpds()
    cpds[0]["dG"] = None                                    # break Vina for the pair
    pairs, fps = find_cliff_pairs(cpds, sim_threshold=0.7, dpic50_threshold=1.0)
    scored = score_cliff_pairs(pairs, cpds, fps)
    assert scored["vina"] == [] and scored["twod"] == []    # paired set stays aligned


def test_bootstrap_paired_deterministic_and_shaped():
    vina = [1, 1, 1, 0, 1, 1, 0, 1, 1, 1]                   # Vina ~0.8 accuracy
    twod = [0, 1, 0, 0, 1, 0, 1, 0, 0, 0]                   # 2D ~0.3 (worse than chance on cliffs)
    a = bootstrap_paired(vina, twod, n=2000, seed=0)
    b = bootstrap_paired(vina, twod, n=2000, seed=0)
    assert a == b                                          # deterministic under a fixed seed
    assert a["n_pairs"] == 10
    assert a["vina_accuracy"]["point"] == 0.8
    assert "ci95" in a["vina_accuracy"] and "paired_vina_minus_2d" in a
    assert isinstance(a["vina_beats_2d"], bool)


def test_bootstrap_empty_is_labelled():
    out = bootstrap_paired([], [])
    assert out["n_pairs"] == 0
