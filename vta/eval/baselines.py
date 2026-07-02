"""Parameter-free enrichment baselines (Phase 11 WI-1).

A structure-based method has demonstrated skill only if it beats trivial baselines on the
SAME actives/decoys split (Wallach & Heifets 2018, J. Chem. Inf. Model. 58:916; Chen 2019;
Sieg 2019). Two baselines:

  * 2D similarity — rank each compound by its maximum ECFP4 (radius 2, 2048-bit) Tanimoto to
    the KNOWN ACTIVES, leave-one-out so an active is never compared to itself. A docking method
    that cannot beat this is merely re-discovering 2D chemical memorization.
  * random — a seeded random ranking; reported as a full bootstrap distribution, not one shuffle.

Both return per-compound scores (higher = more active-like) aligned with the input order, so
they plug straight into `metrics.score_metric` and `significance.paired_bootstrap_delta`.
RDKit only.
"""
from __future__ import annotations

import numpy as np


def _ecfp4(smiles_list: list[str], radius: int = 2, nbits: int = 2048):
    """ECFP4 bit vectors (or None per unparseable SMILES)."""
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nbits)
    fps = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi) if smi else None
        fps.append(gen.GetFingerprint(mol) if mol is not None else None)
    return fps


def baseline_2d_similarity(smiles_list: list[str], labels: list[int],
                           radius: int = 2, nbits: int = 2048) -> list[float]:
    """Max leave-one-out ECFP4 Tanimoto of each compound to the active set.

    For a compound that is itself an active, it is excluded from its own comparison set
    (leave-one-out) so the score reflects similarity to OTHER actives, not to itself.
    """
    from rdkit import DataStructs

    if len(smiles_list) != len(labels):
        raise ValueError("smiles_list and labels must be the same length")
    fps = _ecfp4(smiles_list, radius, nbits)
    active_idx = [i for i, l in enumerate(labels) if l == 1]
    scores = []
    for i, fp in enumerate(fps):
        if fp is None:
            scores.append(0.0)
            continue
        others = [fps[j] for j in active_idx if j != i and fps[j] is not None]
        if not others:
            scores.append(0.0)
            continue
        sims = DataStructs.BulkTanimotoSimilarity(fp, others)
        scores.append(float(max(sims)))
    return scores


def random_metric_distribution(labels: list[int], metric: str, *,
                               n: int = 10000, seed: int = 0, alpha: float = 20.0) -> dict:
    """Null distribution of `metric` under random ranking (WI-1b: full distribution).

    Assigns random scores `n` times and recomputes the metric; returns median + 95% CI. This
    is the chance-level reference the real methods must clear.
    """
    from vta.eval.metrics import score_metric

    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n):
        vals.append(score_metric(rng.random(len(labels)).tolist(), labels, metric, alpha))
    v = np.asarray(vals)
    return {
        "metric": metric,
        "median": round(float(np.median(v)), 4),
        "ci95": [round(float(np.percentile(v, 2.5)), 4),
                 round(float(np.percentile(v, 97.5)), 4)],
        "n_resamples": n,
        "seed": seed,
    }


def bootstrap_metric_ci(scores: list[float], labels: list[int], metric: str, *,
                        n: int = 10000, seed: int = 0, alpha: float = 20.0) -> dict:
    """Bootstrap median + 95% CI of `metric` for one scoring (resample compounds)."""
    from vta.eval.metrics import score_metric

    s, lab = np.asarray(scores, float), np.asarray(labels, int)
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n):
        idx = rng.integers(0, len(lab), size=len(lab))
        ls = lab[idx].tolist()
        if sum(ls) == 0 or sum(ls) == len(ls):
            continue
        vals.append(score_metric(s[idx].tolist(), ls, metric, alpha))
    if not vals:
        return {"metric": metric, "median": None, "ci95": [None, None], "n_effective": 0}
    v = np.asarray(vals)
    return {
        "metric": metric,
        "point": round(score_metric(scores, labels, metric, alpha), 4),
        "median": round(float(np.median(v)), 4),
        "ci95": [round(float(np.percentile(v, 2.5)), 4),
                 round(float(np.percentile(v, 97.5)), 4)],
        "n_effective": len(vals),
        "seed": seed,
    }
