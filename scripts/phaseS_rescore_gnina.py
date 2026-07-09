"""Phase S — GNINA CNN rescoring of activity-cliff poses (pilot: strict Tanimoto ≥ 0.9 cliffs).

Wires a learned rescorer (GNINA's CNNaffinity, McNutt 2021/2025) onto the powered cliff
benchmark. The decisive question: does a CNN rescorer beat BOTH the trivial 2D-kNN QSAR AND
chance on activity cliffs, and does it fix Vina's below-chance *anti-correlation* with potency?

Constraints on this machine (arm64 macOS): no native GNINA; the gnina/gnina Docker image is
amd64 (runs under emulation) and each invocation reloads the CNN ensemble (~60 s). So instead of
one docker call per pose, all pilot poses are combined into a SINGLE multi-ligand SDF and scored
in ONE `--score_only` call — the model loads once. Poses come from the committed Vina docking
(Job A + Phase-S), rescored against the same 7L11 receptor. GNINA CNNaffinity (predicted pK,
higher = more potent) is cached and used as a third method alongside Vina −ΔG and the 2D-kNN QSAR,
scored on the identical strict-cliff pair set through the paired-bootstrap gate.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

IDS_FILE = Path("outputs/phaseS/.strict_cliff_ids.json")
SDF = Path("outputs/phaseS/strict_cliff_poses.sdf")
GNINA_CACHE = Path("outputs/phaseS/.gnina_cache.json")
RECEPTOR = "structures/mpro_30to1_MPRO_primary_rec.pdbqt.receptor.pdb"
REPO = Path(__file__).resolve().parent.parent


def _pose_path(cid: str) -> Optional[str]:
    for pat in (f"structures/phaseS_{cid}_pose.pdbqt",
                f"structures/mpro_30to1_{cid}_pose.pdbqt",
                f"structures/mpro_30to1_{cid}_p1_primary.pdbqt"):
        if os.path.exists(pat):
            return pat
    return None


def build_sdf() -> List[str]:
    """Combine each cliff-member's best Vina pose into one multi-ligand SDF (title = compound id).

    Returns the ordered list of ids actually written (those with a readable pose), so the GNINA
    output — scored in file order — can be mapped back deterministically.
    """
    from openbabel import pybel
    ids = json.loads(IDS_FILE.read_text())["ids"]
    SDF.parent.mkdir(parents=True, exist_ok=True)
    written: List[str] = []
    out = pybel.Outputfile("sdf", str(SDF), overwrite=True)
    for cid in ids:
        p = _pose_path(cid)
        if not p:
            continue
        try:
            mol = next(pybel.readfile("pdbqt", p))   # best (first) model = mode 1
        except Exception:
            continue
        mol.title = cid
        out.write(mol)
        written.append(cid)
    out.close()
    return written


def run_gnina(sdf: Path) -> str:
    """One emulated GNINA --score_only call over the whole SDF (model loads once)."""
    cmd = ["docker", "run", "--rm", "--platform", "linux/amd64",
           "-v", f"{REPO}:/work", "-w", "/work", "gnina/gnina",
           "gnina", "-r", RECEPTOR, "-l", str(sdf),
           "--score_only", "--cnn_scoring", "rescore"]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    (SDF.parent / "gnina_rescore.log").write_text(p.stdout + "\n=== STDERR ===\n" + p.stderr)
    return p.stdout


def parse_cnnaffinity(output: str, ids: List[str]) -> Dict[str, float]:
    """Map per-ligand CNNaffinity (in file order) back to compound ids.

    GNINA prints one block per scored ligand; we read CNNaffinity values in order and zip with
    the ids we wrote. If counts disagree we keep the aligned prefix and label the rest missing.
    """
    vals = [float(m) for m in re.findall(r"CNNaffinity:\s*(-?\d+\.?\d*)", output)]
    out: Dict[str, float] = {}
    for cid, v in zip(ids, vals):
        out[cid] = v
    return out


def _strict_compounds() -> list:
    """Load strict-cliff compound universe with pic50 + dG + gnina CNNaffinity."""
    from vta.eval.cliffs import pic50_from_um
    s = json.loads(Path("vta/data/mpro/mpro_stratified.json").read_text())
    dG: Dict[str, float] = {}
    for p in ("outputs/phase11/.mpro_30to1_dgcache.json", "outputs/phaseS/.cliff_dgcache.json"):
        for k, v in json.loads(Path(p).read_text()).items():
            if isinstance(v, (int, float)):
                dG.setdefault(k, v)
    cnn = json.loads(GNINA_CACHE.read_text()) if GNINA_CACHE.exists() else {}
    cpds = []
    for key in ("non_covalent_actives", "inactives"):
        for r in s.get(key, []):
            pi = pic50_from_um((r.get("measure") or {}).get("ic50_um"))
            smi = r.get("canonical_smiles") or r.get("smiles")
            if pi is not None and smi and r["id"] in dG:
                cpds.append({"id": r["id"], "smiles": smi, "pic50": round(pi, 4),
                             "dG": dG[r["id"]], "cnn": cnn.get(r["id"])})
    return cpds


def score_three_way() -> dict:
    """Score GNINA vs Vina vs 2D-kNN on the strict-cliff pairs where ALL three are available."""
    from vta.eval.cliffs import (bootstrap_paired, find_cliff_pairs, _knn_pic50)
    cpds = _strict_compounds()
    pairs, fps = find_cliff_pairs(cpds, sim_threshold=0.9, dpic50_threshold=1.0)
    vina, gnina, twod = [], [], []
    for (i, j, _s, _d) in pairs:
        ci, cj = cpds[i], cpds[j]
        if ci["dG"] is None or cj["dG"] is None or ci["cnn"] is None or cj["cnn"] is None:
            continue                                  # keep the pair set common to all 3 methods
        mp = i if ci["pic50"] > cj["pic50"] else j
        vina.append(1 if (i if ci["dG"] < cj["dG"] else j) == mp else 0)
        gnina.append(1 if (i if ci["cnn"] > cj["cnn"] else j) == mp else 0)
        pi = _knn_pic50(cpds, fps, i, exclude={i, j}) or -1e9
        pj = _knn_pic50(cpds, fps, j, exclude={i, j}) or -1e9
        twod.append(1 if (i if pi > pj else j) == mp else 0)
    g_vs_2d = bootstrap_paired(gnina, twod, seed=0)
    g_vs_vina = bootstrap_paired(gnina, vina, seed=1)
    demonstrated = bool(g_vs_2d.get("n_pairs") and g_vs_2d["vina_accuracy"]["ci95"][0] > 0.5
                        and g_vs_2d["vina_beats_2d"])   # 'vina_accuracy' key = method A = GNINA
    payload = {
        "phase": "S — GNINA CNN rescoring on strict activity cliffs (Tanimoto ≥ 0.9)",
        "n_pairs_scored": len(gnina),
        "n_poses_rescored": sum(c["cnn"] is not None for c in cpds),
        "gnina_vs_2dknn": g_vs_2d, "gnina_vs_vina": g_vs_vina,
        "accuracies": {
            "gnina": g_vs_2d.get("vina_accuracy"), "twod_knn": g_vs_2d.get("twod_knn_accuracy"),
            "vina": g_vs_vina.get("twod_knn_accuracy")},
        "gnina_beats_2d_and_chance": demonstrated,
        "verdict": ("GNINA beats the 2D-kNN QSAR AND chance on strict cliffs — a learned rescorer "
                    "demonstrates structure-based skill in the fair arena." if demonstrated else
                    "GNINA does NOT clear both the 2D-kNN QSAR and chance on strict cliffs "
                    "(pilot); learned rescoring did not demonstrate structure-based skill here."),
        "caveats": [
            "Pilot: strict Tanimoto≥0.9 cliffs only (small n); 2D-kNN is unusually strong here "
            "(~0.84), the hardest bar. A powered run needs the full 1,193-pair set (GNINA absent "
            "on this arm64 box without emulation).",
            "GNINA CNNaffinity rescored committed Vina poses (--score_only), amd64 under emulation.",
        ],
    }
    Path("outputs/phaseS/gnina_cliff_rescore.json").write_text(json.dumps(payload, indent=2))
    return payload


def rescore() -> Dict[str, float]:
    ids = build_sdf()
    print(f"combined SDF: {len(ids)} poses -> {SDF}", flush=True)
    stdout = run_gnina(SDF)
    cnn = parse_cnnaffinity(stdout, ids)
    GNINA_CACHE.write_text(json.dumps(cnn, indent=2))
    print(f"GNINA CNNaffinity parsed for {len(cnn)}/{len(ids)} poses -> {GNINA_CACHE}", flush=True)
    return cnn


if __name__ == "__main__":
    rescore()
