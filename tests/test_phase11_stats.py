"""Hermetic tests for Phase 11 statistics: score metrics, baselines, paired bootstrap."""
from __future__ import annotations

from vta.eval.baselines import (
    baseline_2d_similarity,
    bootstrap_metric_ci,
    random_metric_distribution,
)
from vta.eval.metrics import score_metric
from vta.eval.significance import holm_bonferroni, paired_bootstrap_delta

PERFECT = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]      # actives on top
INVERTED = list(reversed(PERFECT))              # actives on bottom
LABELS = [1, 1, 1, 1, 1, 0, 0, 0, 0, 0]


def test_score_metric_orientation():
    assert score_metric(PERFECT, LABELS, "roc_auc") == 1.0
    assert score_metric(INVERTED, LABELS, "roc_auc") == 0.0
    assert score_metric(PERFECT, LABELS, "bedroc") > score_metric(INVERTED, LABELS, "bedroc")


def test_paired_bootstrap_detects_real_difference():
    out = paired_bootstrap_delta(PERFECT, INVERTED, LABELS, "roc_auc", n=2000, seed=1)
    assert out["favours_a"] is True
    assert out["ci95"][0] > 0          # CI strictly above 0
    assert out["p_gt0"] > 0.97


def test_paired_bootstrap_identical_methods_not_significant():
    # Same scores → Δ ≡ 0 → not significant, does not favour either method.
    out = paired_bootstrap_delta(PERFECT, PERFECT, LABELS, "bedroc", n=1000, seed=2)
    assert out["significant"] is False
    assert out["favours_a"] is False


def test_paired_bootstrap_is_seed_deterministic():
    a = paired_bootstrap_delta(PERFECT, INVERTED, LABELS, "logauc", n=1000, seed=7)
    b = paired_bootstrap_delta(PERFECT, INVERTED, LABELS, "logauc", n=1000, seed=7)
    assert a["median_delta"] == b["median_delta"] and a["ci95"] == b["ci95"]


def test_random_baseline_is_chance_level():
    out = random_metric_distribution(LABELS, "roc_auc", n=3000, seed=3)
    assert 0.4 < out["median"] < 0.6   # random ranking ≈ chance


def test_2d_similarity_separates_actives_from_distinct_decoy():
    # Two identical actives + one structurally distinct decoy; LOO similarity should rank
    # the actives (each near-identical to the other active) above the decoy.
    smiles = ["CCO", "CCO", "c1ccccc1"]
    labels = [1, 1, 0]
    sims = baseline_2d_similarity(smiles, labels)
    assert sims[0] == 1.0 and sims[1] == 1.0   # leave-one-out: identical other active
    assert sims[2] < sims[0]                    # distinct decoy scores lower
    assert score_metric(sims, labels, "roc_auc") == 1.0


def test_bootstrap_metric_ci_shape():
    out = bootstrap_metric_ci(PERFECT, LABELS, "bedroc", n=1000, seed=4)
    assert out["ci95"][0] <= out["median"] <= out["ci95"][1]
    assert out["n_effective"] > 0


def test_holm_bonferroni_step_down():
    # Smallest p passes the strictest threshold; a large p fails and blocks later ones.
    out = holm_bonferroni({"bedroc": 0.001, "logauc": 0.04, "roc_auc": 0.6, "ef1": 0.5},
                          family_alpha=0.05)
    assert out["bedroc"]["reject"] is True
    assert out["roc_auc"]["reject"] is False
    # step-down: once roc_auc (p=0.6) fails, ef1 cannot be rejected either
    assert out["ef1"]["reject"] is False
