"""Paired-bootstrap significance testing for enrichment method comparisons (Phase 11 WI-3).

The prior pipeline decided significance by asking whether two methods' marginal confidence
intervals overlapped. That test is invalid and over-conservative: two estimates can have
overlapping CIs yet a paired difference that is reliably non-zero (Schenker & Gentleman 2001,
Am. Stat. 55:182; Cumming 2009, Stat. Med. 28:205). Because both methods are scored on the
SAME compounds, the correct test resamples the shared, aligned compound list once per
bootstrap iteration and recomputes BOTH methods' metrics on that same resample, collecting the
paired difference Δ = metric_A − metric_B.

`paired_bootstrap_delta` is the single entry point; `holm_bonferroni` corrects across the
metric family so multiple comparisons do not inflate the false-positive rate.
"""
from __future__ import annotations

import numpy as np

from vta.eval.metrics import score_metric


def paired_bootstrap_delta(scores_a: list[float], scores_b: list[float],
                           labels: list[int], metric: str, *,
                           n: int = 10000, seed: int = 0, alpha: float = 20.0) -> dict:
    """Paired bootstrap of Δ = metric(A) − metric(B) over a shared compound list.

    A and B are two rankings of the SAME compounds (e.g. Vina vs 2D-similarity); `labels`
    are the shared active/decoy labels. Higher score = more active-like. Returns the median Δ,
    its 95% percentile CI, and P(Δ > 0). `significant` is True iff the 95% CI is strictly on
    one side of 0 (Δ reliably non-zero); `favours_a` iff that side is positive.
    """
    n_c = len(labels)
    if not (len(scores_a) == len(scores_b) == n_c):
        raise ValueError("scores_a, scores_b, labels must be the same length")
    sa, sb, lab = np.asarray(scores_a, float), np.asarray(scores_b, float), np.asarray(labels, int)
    rng = np.random.default_rng(seed)

    deltas, a_vals, b_vals = [], [], []
    for _ in range(n):
        idx = rng.integers(0, n_c, size=n_c)
        ls = lab[idx].tolist()
        if sum(ls) == 0 or sum(ls) == len(ls):
            continue  # metric undefined without both classes in the resample
        ma = score_metric(sa[idx].tolist(), ls, metric, alpha)
        mb = score_metric(sb[idx].tolist(), ls, metric, alpha)
        deltas.append(ma - mb)
        a_vals.append(ma)
        b_vals.append(mb)

    if not deltas:
        return {"metric": metric, "n_effective": 0, "median_delta": None,
                "ci95": [None, None], "p_gt0": None, "significant": False, "favours_a": False}
    d = np.asarray(deltas)
    lo, hi = float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))
    p_gt0 = float((d > 0).mean())
    significant = lo > 0 or hi < 0
    return {
        "metric": metric,
        "n_effective": len(deltas),
        "median_delta": round(float(np.median(d)), 4),
        "ci95": [round(lo, 4), round(hi, 4)],
        "p_gt0": round(p_gt0, 4),
        "median_a": round(float(np.median(a_vals)), 4),
        "median_b": round(float(np.median(b_vals)), 4),
        "significant": bool(significant),          # CI strictly excludes 0
        "favours_a": bool(significant and lo > 0),  # A reliably beats B
        "seed": seed,
        "n_resamples": n,
    }


def holm_bonferroni(pvalues: dict[str, float], family_alpha: float = 0.05) -> dict:
    """Holm–Bonferroni step-down correction across a family of comparisons.

    `pvalues` maps a label (e.g. metric name) to a one-sided p-value. Returns, per label, the
    Holm threshold it was tested against and whether it is rejected (significant after
    correction). Controls family-wise error without assuming independence.
    """
    items = sorted(pvalues.items(), key=lambda kv: kv[1])
    m = len(items)
    out, still_rejecting = {}, True
    for rank, (label, p) in enumerate(items):
        threshold = family_alpha / (m - rank)
        reject = still_rejecting and p is not None and p <= threshold
        if not reject:
            still_rejecting = False  # step-down: once one fails, all later fail
        out[label] = {"p": p, "holm_threshold": round(threshold, 5), "reject": bool(reject)}
    return out
