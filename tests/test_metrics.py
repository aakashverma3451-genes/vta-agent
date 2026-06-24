"""Enrichment-metric tests (SPEC #5) — pure, hermetic (no docking, no network)."""
from __future__ import annotations

from vta.eval.metrics import (
    bedroc, enrichment_factor, enrichment_report, roc_auc,
)

# 100 compounds, 10 actives (Ra = 0.1 → EF ceiling = 1/Ra = 10).
N, A = 100, 10
PERFECT = [1] * A + [0] * (N - A)                 # all actives ranked first
WORST = [0] * (N - A) + [1] * A                   # all actives ranked last
SPREAD = [1 if i % 10 == 0 else 0 for i in range(N)]   # one active per decile
SCORES = list(range(N, 0, -1))                    # strictly decreasing → matches order


# ── enrichment factor ────────────────────────────────────────────────────────
def test_ef_perfect_hits_ceiling():
    # Whole top slice is actives → EF = 1/Ra = 10 at every frac up to the active count.
    assert enrichment_factor(PERFECT, 0.01) == 10.0
    assert enrichment_factor(PERFECT, 0.05) == 10.0
    assert enrichment_factor(PERFECT, 0.10) == 10.0


def test_ef_worst_is_zero():
    assert enrichment_factor(WORST, 0.01) == 0.0
    assert enrichment_factor(WORST, 0.10) == 0.0


def test_ef_spread_is_about_one():
    # Evenly spread actives → top 10% holds exactly its fair share → EF ≈ 1.
    assert abs(enrichment_factor(SPREAD, 0.10) - 1.0) < 1e-6


def test_ef_degenerate_inputs_are_zero():
    assert enrichment_factor([], 0.1) == 0.0
    assert enrichment_factor([0, 0, 0], 0.1) == 0.0       # no actives
    assert enrichment_factor([1, 1], 0.0) == 0.0          # non-positive frac


# ── ROC-AUC ──────────────────────────────────────────────────────────────────
def test_auc_perfect_is_one():
    assert roc_auc(SCORES, PERFECT) == 1.0


def test_auc_inverted_is_zero():
    assert roc_auc(SCORES, WORST) == 0.0


def test_auc_random_is_about_half():
    assert abs(roc_auc(SCORES, SPREAD) - 0.5) < 0.1


def test_auc_all_ties_is_half():
    # Identical scores → every comparison is a tie → 0.5.
    assert roc_auc([1.0] * 4, [1, 0, 1, 0]) == 0.5


def test_auc_undefined_when_one_class_empty():
    assert roc_auc([3, 2, 1], [1, 1, 1]) == 0.0          # no decoys


# ── BEDROC ───────────────────────────────────────────────────────────────────
def test_bedroc_perfect_is_one():
    assert bedroc(PERFECT) > 0.99


def test_bedroc_worst_is_zero():
    assert bedroc(WORST) < 0.01


def test_bedroc_spread_is_low():
    # A spread/random ranking has no early recognition → small BEDROC (near Ra).
    assert bedroc(SPREAD) < 0.5


def test_bedroc_in_unit_interval():
    for labels in (PERFECT, WORST, SPREAD):
        assert 0.0 <= bedroc(labels) <= 1.0


def test_bedroc_rewards_earlier_recognition():
    # Same number of actives, but earlier placement must score higher.
    early = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0]
    late = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1]
    assert bedroc(early) > bedroc(late)


# ── enrichment_report (the benchmark core) ────────────────────────────────────
def _entries(actives_scores, decoy_scores):
    e = [{"name": f"a{i}", "score": s, "positive_control": True}
         for i, s in enumerate(actives_scores)]
    e += [{"name": f"d{i}", "score": s, "positive_control": False}
          for i, s in enumerate(decoy_scores)]
    return e


def test_report_perfect_separation():
    rep = enrichment_report(_entries([0.9, 0.8, 0.7], [0.3, 0.2, 0.1]))
    assert rep["n"] == 6 and rep["n_actives"] == 3 and rep["n_decoys"] == 3
    assert rep["roc_auc"] == 1.0
    assert rep["bedroc"] > 0.9
    assert rep["ef"]["EF20%"] == 2.0          # ceiling 1/Ra = 1/0.5
    assert [r["label"] for r in rep["ranking"]] == ["active"] * 3 + ["decoy"] * 3


def test_report_sorts_and_labels_by_score():
    # Unordered input must be ranked best-first with correct labels.
    rep = enrichment_report(_entries([0.4], [0.9, 0.1]))
    assert [r["name"] for r in rep["ranking"]] == ["d0", "a0", "d1"]
    assert rep["roc_auc"] == 0.5              # one active between two decoys


def test_report_inverted_is_zero_auc():
    rep = enrichment_report(_entries([0.1, 0.2], [0.8, 0.9]))
    assert rep["roc_auc"] == 0.0
    assert rep["bedroc"] < 0.1


# ── decoy seam (committed cache) ──────────────────────────────────────────────
def test_fetch_decoys_returns_committed_pool():
    from vta.data import decoys
    out = decoys.fetch_decoys(["CCO"], n=50)
    assert out and all(isinstance(s, str) for s in out)


def test_fetch_decoys_excludes_actives_and_caps():
    from vta.data import decoys
    pool = [d["smiles"] for d in decoys.load_decoy_pool()]
    assert pool, "committed demo decoy pool should be non-empty"
    # An active equal to a pool member is excluded; n caps the per-active count.
    out = decoys.fetch_decoys([pool[0]], n=2)
    assert pool[0] not in out
    assert len(out) <= 2


def test_fetch_decoys_none_when_pool_absent(monkeypatch, tmp_path):
    from vta.data import decoys
    monkeypatch.setattr(decoys, "_DECOY_DIR", str(tmp_path / "nope"))
    assert decoys.fetch_decoys(["CCO"]) is None
