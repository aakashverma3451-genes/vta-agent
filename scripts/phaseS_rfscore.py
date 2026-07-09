"""Phase S — in-house RF-Score learned rescorer on the powered activity-cliff benchmark.

Asks whether a LEARNED scorer (RF-Score, Ballester & Mitchell 2010) beats the trivial 2D-kNN
QSAR on Mpro activity cliffs — and whether it fixes Vina's below-chance anti-correlation. Built
on the working stack (openbabel + numpy + scikit-learn) because GNINA/RTMScore/ODDT are all
un-installable on this arm64 box.

Pipeline: RF-Score contact-count features from every committed Vina pose → RandomForest with
**scaffold-clustered leave-out CV** (no analog leakage) → out-of-fold pIC50 predictions →
scored on the SAME 1,193 cliff pairs as Vina and the 2D-kNN QSAR, through the paired-bootstrap
gate. This is a TARGET-SPECIFIC model (trained on Moonshot Mpro), so a win is a target-specific
result, not a general-transfer claim — stated plainly in the verdict.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from vta.eval.cliffs import (bootstrap_paired, find_cliff_pairs, pic50_from_um,
                             rank_by_predictor, score_cliff_pairs)
from vta.eval.rfscore import murcko_scaffold, read_atoms, rfscore_features, scaffold_cv_predict

RECEPTOR_PDB = "structures/mpro_30to1_MPRO_primary_rec.pdbqt.receptor.pdb"
STRAT = Path("vta/data/mpro/mpro_stratified.json")
OUT = Path("outputs/phaseS/rfscore_cliff.json")
import os


def _pose_path(cid: str) -> Optional[str]:
    for pat in (f"structures/phaseS_{cid}_pose.pdbqt",
                f"structures/mpro_30to1_{cid}_pose.pdbqt",
                f"structures/mpro_30to1_MPRO_{cid}_p1_primary.pdbqt"):
        if os.path.exists(pat):
            return pat
    return None


def _load_docked() -> List[dict]:
    # Single source of truth: reuse the committed cliff benchmark's universe so the 2D-kNN
    # baseline and Vina are provably identical to phaseS_activity_cliffs; just attach pose + rf.
    from scripts.phaseS_activity_cliffs import _load_compounds as _cliff_universe
    cpds = _cliff_universe()
    for c in cpds:
        c["pose"] = _pose_path(c["id"])
        c["rf"] = None
    return cpds


def run() -> dict:
    receptor = read_atoms(RECEPTOR_PDB, "pdb")
    cpds = _load_docked()

    # RF-Score features from each pose
    feated = []
    for c in cpds:
        if c["pose"]:
            f = rfscore_features(receptor, c["pose"])
            if f is not None:
                c["feat"] = f
                feated.append(c)
    X = np.asarray([c["feat"] for c in feated])
    y = np.asarray([c["pic50"] for c in feated])
    groups = [murcko_scaffold(c["smiles"]) or f"__id__{c['id']}" for c in feated]
    n_scaffolds = len(set(groups))
    preds = scaffold_cv_predict(X, y, groups, n_splits=5, seed=0)
    for c, p in zip(feated, preds):
        c["rf"] = None if np.isnan(p) else float(p)

    # Score Vina + 2D-kNN via the SAME code path as the committed cliff benchmark (so the 2D
    # baseline is provably identical), then add RF on that exact pair set via rank_by_predictor.
    pairs, fps = find_cliff_pairs(cpds, sim_threshold=0.7, dpic50_threshold=1.0)
    scored = score_cliff_pairs(pairs, cpds, fps)                 # {"vina","twod","pairs"}
    common = set(scored["pairs"])
    rf_flags, rf_pairs = rank_by_predictor(pairs, cpds, lambda c: c.get("rf"), restrict=common)
    # align Vina/2D to exactly the pairs RF could score (all, unless a pose was unreadable)
    keep = set(rf_pairs)
    idx = [n for n, pr in enumerate(scored["pairs"]) if pr in keep]
    vina = [scored["vina"][n] for n in idx]
    twod = [scored["twod"][n] for n in idx]
    rf = rf_flags

    rf_vs_2d = bootstrap_paired(rf, twod, seed=0)      # method A = RF-Score, B = 2D-kNN
    rf_vs_vina = bootstrap_paired(rf, vina, seed=1)
    demonstrated = bool(rf_vs_2d.get("n_pairs")
                        and rf_vs_2d["vina_accuracy"]["ci95"][0] > 0.5   # RF acc CI excludes 0.5
                        and rf_vs_2d["vina_beats_2d"])                    # RF beats 2D (paired CI>0)

    payload = {
        "phase": "S — in-house RF-Score learned rescorer on powered activity cliffs",
        "method": "RF-Score v1 contact-count features + RandomForest; scaffold-clustered leave-out CV",
        "n_compounds_featurized": len(feated),
        "n_scaffold_groups": n_scaffolds,
        "n_cliff_pairs_scored": len(rf),
        "accuracies": {
            "rfscore": rf_vs_2d.get("vina_accuracy"),
            "twod_knn": rf_vs_2d.get("twod_knn_accuracy"),
            "vina": rf_vs_vina.get("twod_knn_accuracy"),
            "chance": 0.5,
        },
        "rfscore_vs_2dknn": rf_vs_2d, "rfscore_vs_vina": rf_vs_vina,
        "rfscore_beats_2d_and_chance": demonstrated,
        "verdict": (
            "TARGET-SPECIFIC learned rescoring beats the 2D-kNN QSAR AND chance on activity cliffs "
            "(scaffold-CV, out-of-fold) — a learned scorer shows cliff-ranking skill here." if demonstrated
            else "The learned RF-Score rescorer does NOT beat both the 2D-kNN QSAR and chance on "
                 "activity cliffs (scaffold-CV, out-of-fold) — learned rescoring did not demonstrate "
                 "cliff-ranking skill over the trivial baseline."),
        "caveats": [
            "TARGET-SPECIFIC: RF-Score is trained on Moonshot Mpro itself (out-of-fold under "
            "scaffold-clustered CV — no analog leakage), so this is NOT a general-transfer claim. "
            "A PDBbind-pretrained scorer (GNINA/RTMScore) on a native box answers that separately.",
            "RF-Score v1 (coarse contact counts) is weaker than modern CNN/GNN scorers; a null "
            "here is suggestive, not proof that GNINA would also fail on these cliffs.",
            "Poses are the committed Vina docks; the same 1,193-pair set and 2D-kNN baseline as "
            "the powered cliff benchmark, through the identical paired gate.",
        ],
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2))
    _md(payload)
    return payload


def _md(p: dict) -> None:
    a = p["accuracies"]
    def acc(x):
        return f"{x['median']} {x['ci95']}" if x else "—"
    lines = [f"# Phase S — in-house RF-Score learned rescorer on activity cliffs", "",
             p["method"], "",
             f"- compounds featurized: {p['n_compounds_featurized']} across "
             f"{p['n_scaffold_groups']} scaffold groups",
             f"- cliff pairs scored (RF+Vina+2D all available): {p['n_cliff_pairs_scored']}", "",
             "## Cliff-pair ranking accuracy (median [95% CI]); chance = 0.50", "",
             "| Method | accuracy |", "|---|---|",
             f"| **RF-Score (learned)** | {acc(a['rfscore'])} |",
             f"| 2D-kNN QSAR | {acc(a['twod_knn'])} |",
             f"| Vina (−ΔG) | {acc(a['vina'])} |", "",
             f"Paired RF−2D: Δ {p['rfscore_vs_2dknn']['paired_vina_minus_2d']['median']} "
             f"{p['rfscore_vs_2dknn']['paired_vina_minus_2d']['ci95']}; "
             f"RF−Vina: Δ {p['rfscore_vs_vina']['paired_vina_minus_2d']['median']} "
             f"{p['rfscore_vs_vina']['paired_vina_minus_2d']['ci95']}.", "",
             f"**{p['verdict']}**", "", "## Caveats"] + [f"- {c}" for c in p["caveats"]] + [""]
    OUT.with_suffix(".md").write_text("\n".join(lines))


def main() -> None:
    p = run()
    a = p["accuracies"]
    print(f"featurized={p['n_compounds_featurized']} scaffolds={p['n_scaffold_groups']} "
          f"pairs={p['n_cliff_pairs_scored']}")
    print(f"RF {a['rfscore']['median']} {a['rfscore']['ci95']} | 2D {a['twod_knn']['median']} "
          f"{a['twod_knn']['ci95']} | Vina {a['vina']['median']} {a['vina']['ci95']}")
    print("RF beats 2D+chance:", p["rfscore_beats_2d_and_chance"])
    print(p["verdict"])


if __name__ == "__main__":
    main()
