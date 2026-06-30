"""vta.eval.metrics — virtual-screening enrichment metrics (SPEC #5), pure & hermetic.

A 9-ligand control recovery check ("do the known inhibitors land top-5?") is a sanity
test, not a benchmark. To quantify screening power you dock a curated ACTIVE set plus
property-matched DECOYS and measure how strongly the score separates them. This module
provides the three metrics a reviewer expects, with NO heavy dependencies so they stay
fully unit-testable:

    enrichment_factor(labels, frac)  — EF at the top `frac` (e.g. 0.01 = EF1%): how many
        more actives you find in the top frac than random selection would. Lead metric.
    bedroc(labels, alpha)            — Boltzmann-Enhanced Discrimination of ROC (Truchon
        & Bayly 2007, doi:10.1021/ci600426e). Early-recognition metric in [0,1]; alpha=20
        weights the top ~8% (the slice virtual screening actually acts on).
    roc_auc(scores, labels)          — global ranking quality in [0,1]. Reported for
        completeness, but it over-weights the late (VS-irrelevant) tail, so it is NOT the
        headline number for early recognition.

Convention: higher score = more active-like. `labels` for EF/BEDROC are ordered BEST
score first (1 = active, 0 = decoy); `roc_auc` takes raw (scores, labels) and ranks
internally.
"""
from __future__ import annotations

import math
import random


def enrichment_factor(labels: list[int], frac: float) -> float:
    """EF at the top `frac` of a best-first-ranked label list.

    EF = (actives among the top `frac` / size of that top slice)
       / (total actives / total compounds)

    1.0 = no better than random; the ceiling is 1/Ra (Ra = active fraction) when the
    whole top slice is actives. Returns 0.0 for degenerate inputs.
    """
    n = len(labels)
    n_act = sum(labels)
    if n == 0 or n_act == 0 or frac <= 0.0:
        return 0.0
    n_top = max(1, int(round(frac * n)))
    hits_top = sum(labels[:n_top])
    return round((hits_top / n_top) / (n_act / n), 4)


def roc_auc(scores: list[float], labels: list[int]) -> float:
    """ROC-AUC = P(score(active) > score(decoy)), ties counted as 0.5.

    Computed via the Mann–Whitney U identity (no sorting subtleties, ties handled). 1.0
    = every active outranks every decoy; 0.5 = random; 0.0 = perfectly inverted. Returns
    0.0 when one class is empty (AUC undefined).
    """
    pos = [s for s, l in zip(scores, labels) if l == 1]
    neg = [s for s, l in zip(scores, labels) if l == 0]
    if not pos or not neg:
        return 0.0
    wins = 0.0
    for p in pos:
        for q in neg:
            if p > q:
                wins += 1.0
            elif p == q:
                wins += 0.5
    return round(wins / (len(pos) * len(neg)), 4)


def log_auc(labels: list[int], min_fp_rate: float = 0.001) -> float:
    """Log-scaled early AUC for a best-first label list.

    This approximates the common virtual-screening logAUC by integrating TPR over a
    logarithmic FPR axis. It is intentionally dependency-free for hermetic tests.
    """
    n = len(labels)
    n_act = sum(labels)
    n_dec = n - n_act
    if n == 0 or n_act == 0 or n_dec == 0:
        return 0.0
    fp = tp = 0
    points = [(min_fp_rate, 0.0)]
    for label in labels:
        if label:
            tp += 1
        else:
            fp += 1
        fpr = max(min_fp_rate, fp / n_dec)
        tpr = tp / n_act
        points.append((fpr, tpr))
    area = 0.0
    denom = abs(math.log10(min_fp_rate))
    for (x1, y1), (x2, y2) in zip(points, points[1:]):
        lx1, lx2 = math.log10(x1), math.log10(x2)
        area += ((y1 + y2) / 2.0) * (lx2 - lx1)
    return round(max(0.0, min(1.0, area / denom)), 4)


def bedroc(labels: list[int], alpha: float = 20.0) -> float:
    """BEDROC in [0,1] for a best-first-ranked label list (Truchon & Bayly 2007).

    Exponentially weights early ranks (alpha=20 → the top ~8%), so it rewards finding
    actives EARLY — the property virtual screening cares about and ROC-AUC misses.
    ~1.0 = actives concentrated at the top; ~0.0 = actives at the bottom; a random
    ranking gives a low value near Ra. Returns 0.0 when there are no actives or no
    decoys (BEDROC undefined).
    """
    n = len(labels)
    n_act = sum(labels)
    if n == 0 or n_act == 0 or n_act == n:
        return 0.0
    ra = n_act / n                                   # fraction of actives
    # 1-based ranks of the actives in the best-first list.
    ranks = [i + 1 for i, l in enumerate(labels) if l == 1]
    s = sum(math.exp(-alpha * r / n) for r in ranks)
    # Random-expectation normaliser (denominator of the RIE).
    rie_rand = (1.0 - math.exp(-alpha)) / (math.exp(alpha / n) - 1.0)
    rie = s / (ra * rie_rand)
    # Closed-form map RIE → BEDROC ∈ [0,1] (Truchon & Bayly eq. 36).
    fac = ra * math.sinh(alpha / 2.0) / (
        math.cosh(alpha / 2.0) - math.cosh(alpha / 2.0 - alpha * ra))
    offset = 1.0 / (1.0 - math.exp(alpha * (1.0 - ra)))
    return round(rie * fac + offset, 4)


def enrichment_report(entries: list[dict],
                      fractions=(0.01, 0.05, 0.10, 0.20),
                      alpha: float = 20.0) -> dict:
    """Full enrichment summary for scored compounds — pure, the benchmark's core.

    `entries`: [{name, score, positive_control: bool}] (any order). Ranks best-first by
    score and returns EF at each fraction, BEDROC, ROC-AUC, and the ranked table. The
    driver wraps this with I/O + the decoy-bias caveat; tests exercise it directly.
    """
    ranked = sorted(entries, key=lambda e: e["score"], reverse=True)
    labels = [1 if e.get("positive_control") else 0 for e in ranked]
    scores = [e["score"] for e in ranked]
    return {
        "n": len(ranked),
        "n_actives": sum(labels),
        "n_decoys": len(labels) - sum(labels),
        "roc_auc": roc_auc(scores, labels),
        "bedroc": bedroc(labels, alpha),
        "log_auc": log_auc(labels),
        "ef": {f"EF{int(f * 100)}%": enrichment_factor(labels, f) for f in fractions},
        "ranking": [
            {"rank": i + 1, "name": e["name"], "score": e["score"],
             "label": "active" if labels[i] else "decoy"}
            for i, e in enumerate(ranked)
        ],
    }


def _percentile(values: list[float], pct: float) -> float:
    if not values:
        return 0.0
    vals = sorted(values)
    k = (len(vals) - 1) * pct
    lo, hi = math.floor(k), math.ceil(k)
    if lo == hi:
        return vals[int(k)]
    return vals[lo] * (hi - k) + vals[hi] * (k - lo)


def bootstrap_enrichment_report(entries: list[dict], *,
                                n_resamples: int = 1000,
                                seed: int = 42,
                                alpha: float = 20.0) -> dict:
    """Bootstrap CIs for enrichment metrics.

    Resamples compounds with replacement. Degenerate samples lacking either class are
    skipped, because ROC/BEDROC are undefined there. Returns point estimates plus
    median and 95% CI for each primary metric.
    """
    point = enrichment_report(entries, alpha=alpha)
    if not entries:
        return {"point": point, "bootstrap": {}, "n_resamples": 0, "skipped": n_resamples}

    rng = random.Random(seed)
    metrics = {"bedroc": [], "log_auc": [], "roc_auc": [], "EF1%": []}
    skipped = 0
    for _ in range(n_resamples):
        sample = [entries[rng.randrange(len(entries))] for _ in entries]
        labels = [1 if e.get("positive_control") else 0 for e in sample]
        if sum(labels) == 0 or sum(labels) == len(labels):
            skipped += 1
            continue
        rep = enrichment_report(sample, fractions=(0.01,), alpha=alpha)
        metrics["bedroc"].append(rep["bedroc"])
        metrics["log_auc"].append(rep["log_auc"])
        metrics["roc_auc"].append(rep["roc_auc"])
        metrics["EF1%"].append(rep["ef"]["EF1%"])

    boot = {}
    for name, vals in metrics.items():
        boot[name] = {
            "median": round(_percentile(vals, 0.50), 4),
            "ci95": [
                round(_percentile(vals, 0.025), 4),
                round(_percentile(vals, 0.975), 4),
            ],
        }
    return {
        "point": point,
        "bootstrap": boot,
        "n_resamples": n_resamples,
        "n_effective": len(metrics["bedroc"]),
        "skipped": skipped,
        "seed": seed,
    }
