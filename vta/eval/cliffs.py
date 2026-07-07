"""Activity-cliff / matched-molecular-pair scoring (Phase S).

On analog-clustered sets (e.g. COVID Moonshot) 2D-similarity is near-optimal *by design*, so a
docking method losing to it there is expected, not informative (Wallach 2018; the campaign was
an analog-progression effort — Boby 2023). **Activity cliffs invert that**: pairs that are
near-identical in 2D (ECFP4) but differ sharply in potency (ΔpIC50) are the one regime where a
ligand-only 2D method is forced to ~chance, so a structure-based ranker that reads the actual
pocket interaction *can* win — and a win there is meaningful (van Tilborg 2022, MoleculeACE).

This module is pure/RDKit-only and hermetic:
  - `find_cliff_pairs` mines cliff pairs (Tanimoto ≥ sim_threshold AND |ΔpIC50| ≥ dpic50_threshold),
  - `score_cliff_pairs` asks, per pair, whether each method ranks the MORE-POTENT member higher,
  - the 2D baseline is a kNN-pIC50 predictor that **excludes both pair members** (its fairest
    possible estimate — not rigged to lose),
  - `bootstrap_paired` gives median accuracies + 95% CIs and a paired Vina−2D difference, so the
    result plugs into the same honesty gate as every other benchmark.
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from vta.eval.baselines import _ecfp4


def pic50_from_um(ic50_um: Optional[float]) -> Optional[float]:
    """µM IC50 → pIC50 (= -log10[M]). None/non-positive → None (labelled skip, not fabricated)."""
    if ic50_um is None:
        return None
    try:
        v = float(ic50_um)
    except (TypeError, ValueError):
        return None
    return None if v <= 0 else 6.0 - math.log10(v)


def _tanimoto(a, b) -> float:
    from rdkit import DataStructs
    return DataStructs.TanimotoSimilarity(a, b)


def find_cliff_pairs(compounds: List[Dict[str, Any]], *, sim_threshold: float = 0.7,
                     dpic50_threshold: float = 1.0) -> Tuple[List[Tuple[int, int, float, float]], list]:
    """Return [(i, j, tanimoto, |ΔpIC50|)] for cliff pairs, plus the ECFP4 fingerprints.

    A cliff pair: ECFP4 Tanimoto ≥ sim_threshold (2D-similar) AND |ΔpIC50| ≥ dpic50_threshold
    (large potency gap). Compounds must carry `smiles` and a numeric `pic50`.
    """
    fps = _ecfp4([c.get("smiles") or "" for c in compounds])
    n = len(compounds)
    pairs: List[Tuple[int, int, float, float]] = []
    for i in range(n):
        if fps[i] is None or compounds[i].get("pic50") is None:
            continue
        for j in range(i + 1, n):
            if fps[j] is None or compounds[j].get("pic50") is None:
                continue
            sim = _tanimoto(fps[i], fps[j])
            if sim < sim_threshold:
                continue
            dp = abs(compounds[i]["pic50"] - compounds[j]["pic50"])
            if dp >= dpic50_threshold:
                pairs.append((i, j, round(sim, 4), round(dp, 4)))
    return pairs, fps


def _knn_pic50(compounds, fps, target: int, exclude: set, k: int = 3) -> Optional[float]:
    """Predict a compound's pIC50 from its k nearest ECFP4 neighbours (excluding `exclude`)."""
    sims = []
    for m in range(len(compounds)):
        if m in exclude or fps[m] is None or compounds[m].get("pic50") is None:
            continue
        sims.append((_tanimoto(fps[target], fps[m]), compounds[m]["pic50"]))
    if not sims:
        return None
    sims.sort(key=lambda x: x[0], reverse=True)
    top = sims[:k]
    return sum(p for _, p in top) / len(top)


def score_cliff_pairs(pairs, compounds, fps, *, k: int = 3) -> Dict[str, List[int]]:
    """Per pair: does Vina (−ΔG) and the 2D-kNN baseline rank the MORE-POTENT member higher?

    Returns {"vina": [0/1...], "twod": [0/1...]} aligned with `pairs`. A pair needs a numeric
    `dG` on both members to score Vina; pairs missing ΔG are dropped from BOTH lists so the
    comparison stays paired on the identical pair set.
    """
    vina, twod, kept = [], [], []
    for (i, j, _sim, _dp) in pairs:
        di, dj = compounds[i].get("dG"), compounds[j].get("dG")
        if di is None or dj is None:
            continue
        more_potent = i if compounds[i]["pic50"] > compounds[j]["pic50"] else j
        # Vina: predict the more potent as the one with the LOWER (more negative) ΔG.
        vina_pick = i if di < dj else j
        vina.append(1 if vina_pick == more_potent else 0)
        # 2D kNN excluding BOTH pair members (the baseline's fairest estimate).
        pi = _knn_pic50(compounds, fps, i, exclude={i, j}, k=k)
        pj = _knn_pic50(compounds, fps, j, exclude={i, j}, k=k)
        pi = -math.inf if pi is None else pi
        pj = -math.inf if pj is None else pj
        twod_pick = i if pi > pj else j
        twod.append(1 if twod_pick == more_potent else 0)
        kept.append((i, j))
    return {"vina": vina, "twod": twod, "pairs": kept}


def _ci(vals: np.ndarray) -> List[float]:
    return [round(float(np.percentile(vals, 2.5)), 4), round(float(np.percentile(vals, 97.5)), 4)]


def bootstrap_paired(vina: List[int], twod: List[int], *, n: int = 10000, seed: int = 0) -> Dict[str, Any]:
    """Pair-level bootstrap: median accuracy + 95% CI for each method and the paired Vina−2D Δ.

    `significant_vs_chance` = Vina accuracy CI excludes 0.5. `beats_2d` = paired Δ CI excludes 0
    (and favours Vina). Both must hold to claim demonstrated structure-based skill on cliffs.
    """
    va, td = np.array(vina, dtype=float), np.array(twod, dtype=float)
    m = len(va)
    if m == 0:
        return {"n_pairs": 0, "note": "no scorable cliff pairs"}
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, m, size=(n, m))
    va_boot = va[idx].mean(axis=1)
    td_boot = td[idx].mean(axis=1)
    diff = va_boot - td_boot
    va_ci, td_ci, d_ci = _ci(va_boot), _ci(td_boot), _ci(diff)
    return {
        "n_pairs": m,
        "vina_accuracy": {"median": round(float(np.median(va_boot)), 4), "ci95": va_ci,
                          "point": round(float(va.mean()), 4)},
        "twod_knn_accuracy": {"median": round(float(np.median(td_boot)), 4), "ci95": td_ci,
                              "point": round(float(td.mean()), 4)},
        "paired_vina_minus_2d": {"median": round(float(np.median(diff)), 4), "ci95": d_ci,
                                 "p_gt0": round(float((diff > 0).mean()), 4)},
        "chance": 0.5,
        "vina_significant_vs_chance": bool(va_ci[0] > 0.5 or va_ci[1] < 0.5),
        "vina_beats_2d": bool(d_ci[0] > 0.0),
    }
