"""Phase 11 WI-4 Job A — Mpro non-covalent enrichment at a realistic 30:1 ratio.

The frozen headline benchmark is 1:1 (EF1% ceiling 2.0). This re-docks 30 non-covalent
Moonshot actives against 900 experimentally-measured Moonshot inactives (30:1, all real, no
fabrication, no DeepCoy) with the production Vina protocol, then reports enrichment with
bootstrap CIs AND the WI-1 trivial baselines + WI-3 paired test at the new ratio. Writes to a
SEPARATE artifact — the frozen 1:1 headline and locked_benchmark.json are not touched.

Question answered: does Vina's (lack of) enrichment vs 2D-similarity survive a non-degenerate
inactive:active ratio, and does EF become informative once its ceiling is lifted?
"""
from __future__ import annotations

import json
from pathlib import Path

import glob
import re
import subprocess

from vta.eval.baselines import baseline_2d_similarity, bootstrap_metric_ci, random_metric_distribution
from vta.eval.metrics import METRIC_NAMES, bootstrap_enrichment_report, enrichment_report
from vta.eval.significance import holm_bonferroni, paired_bootstrap_delta
from vta.nodes.docking import _prep_ligand, _prep_receptor, _vina_bin
from vta.nodes.pockets import pockets_node
from vta.nodes.structure import structure_node
from vta.state import new_state
from scripts.run_mpro_benchmark import select_library, to_ligand

OUT = Path("outputs/phase11/mpro_30to1_benchmark.json")
DOCK_CACHE = Path("outputs/phase11/.mpro_30to1_dgcache.json")
DOCK_TIMEOUT_S = 180      # hard per-dock cap: a hung Vina search is skipped, not waited on
_LIG_CACHE = "structures/.ligand_pdbqt"
_BOX = 22.0
# Same-assay (Moonshot) inactives cap at 745, so to meet ≥30 inactives/active cleanly we use
# 24 actives × 745 real Moonshot inactives = 31:1 (all one fluorescence assay; no fabrication).
N_ACTIVES, N_INACTIVES, SEED, N_BOOT = 24, 745, 17, 10000


def _vina_dg(text: str) -> float | None:
    for line in text.splitlines():
        m = re.match(r"\s*1\s+(-?\d+\.\d+)", line)
        if m:
            return float(m.group(1))
    return None


def _dock_one(receptor: str, lig_pdbqt: str, center: list, out: str) -> float | None:
    """One Vina dock with a hard timeout; raises subprocess.TimeoutExpired if it hangs."""
    vina = _vina_bin()
    p = subprocess.run(
        [vina, "--receptor", receptor, "--ligand", lig_pdbqt,
         "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
         "--size_x", str(_BOX), "--size_y", str(_BOX), "--size_z", str(_BOX),
         "--exhaustiveness", "8", "--num_modes", "5", "--seed", "42", "--out", out],
        capture_output=True, text=True, timeout=DOCK_TIMEOUT_S)
    return _vina_dg(p.stdout)


def _load_cache() -> dict:
    try:
        return json.loads(DOCK_CACHE.read_text())
    except Exception:
        return {}


def _save_cache(cache: dict) -> None:
    DOCK_CACHE.parent.mkdir(parents=True, exist_ok=True)
    DOCK_CACHE.write_text(json.dumps(cache))


def _seed_cache_from_existing_poses(cache: dict) -> dict:
    """Recover ΔG from poses docked by earlier (interrupted) runs so we resume, not restart."""
    for pose in glob.glob("structures/mpro_30to1_MPRO_*_p1_primary.pdbqt"):
        cid = Path(pose).name[len("mpro_30to1_MPRO_"):-len("_p1_primary.pdbqt")]
        if cid in cache:
            continue
        try:
            for line in Path(pose).read_text().splitlines():
                m = re.search(r"REMARK VINA RESULT:\s*(-?\d+\.\d+)", line)
                if m:
                    cache[cid] = float(m.group(1))
                    break
        except Exception:
            pass
    return cache


def run() -> dict:
    strat = json.loads(Path("vta/data/mpro/mpro_stratified.json").read_text())
    acts, inacts = select_library(strat, N_ACTIVES, N_INACTIVES, SEED)
    library = ([to_ligand(r, True, i) for i, r in enumerate(acts)]
               + [to_ligand(r, False, i + 100000) for i, r in enumerate(inacts)])

    st = new_state("mpro_30to1", "mpro_30to1")
    st["extracted_proteins"] = {"MPRO": {"sequence": "M" * 306, "length_aa": 306, "plddt": None}}
    st = structure_node(st)
    st = pockets_node(st)
    struct = st["structures"]["MPRO"]
    center = st["pockets"]["MPRO"][0]["center"]
    receptor = _prep_receptor(struct["pdb_path"], "structures/mpro_30to1_MPRO_primary_rec")

    # Resumable, hang-capped dock loop. Docking directly (not via docking_node) so we can
    # (a) skip already-docked ligands via a persisted ΔG cache, and (b) hard-cap each dock at
    # DOCK_TIMEOUT_S — a few large/flexible inactives make Vina's search run for tens of
    # minutes and stall the whole sequential screen. A timed-out ligand is skipped (ΔG=None,
    # labelled), never allowed to hang the run.
    cache = _seed_cache_from_existing_poses(_load_cache())
    entries, n_timeout = [], 0
    for i, lig in enumerate(library):
        cid = lig["dock_cache_id"]
        if cid not in cache:
            lp = _prep_ligand(lig, _LIG_CACHE)
            try:
                cache[cid] = (_dock_one(receptor, lp, center,
                                        f"structures/mpro_30to1_{cid}_pose.pdbqt") if lp else None)
            except subprocess.TimeoutExpired:
                cache[cid] = "TIMEOUT"
            _save_cache(cache)
        dG = cache[cid]
        if dG == "TIMEOUT":
            n_timeout += 1
            continue
        if dG is not None:
            entries.append({"name": lig["name"], "score": round(-float(dG), 4),
                            "positive_control": lig["positive_control"],
                            "smiles": lig["smiles"], "dG": dG})
    entries.sort(key=lambda e: e["score"], reverse=True)
    real_vina = bool(entries)

    labels = [1 if e["positive_control"] else 0 for e in entries]
    smiles = [e.get("smiles") or "" for e in entries]
    vina_scores = [e["score"] for e in entries]         # −ΔG
    sim_scores = baseline_2d_similarity(smiles, labels)

    report = enrichment_report(entries)
    bootstrap = bootstrap_enrichment_report(entries, n_resamples=N_BOOT, seed=42)
    baselines = {
        "random": {m: random_metric_distribution(labels, m, n=N_BOOT, seed=11) for m in METRIC_NAMES},
        "2d_sim": {m: bootstrap_metric_ci(sim_scores, labels, m, n=N_BOOT, seed=11) for m in METRIC_NAMES},
        "vina": {m: bootstrap_metric_ci(vina_scores, labels, m, n=N_BOOT, seed=11) for m in METRIC_NAMES},
    }
    vina_vs_2d = {m: paired_bootstrap_delta(vina_scores, sim_scores, labels, m, n=N_BOOT, seed=21)
                  for m in METRIC_NAMES}
    holm = holm_bonferroni({m: (1.0 - vina_vs_2d[m]["p_gt0"]) if vina_vs_2d[m]["p_gt0"] is not None
                            else None for m in METRIC_NAMES})
    structure_skill = bool(vina_vs_2d["bedroc"]["favours_a"] and holm["bedroc"]["reject"])

    n_act = sum(labels)
    payload = {
        "phase": "11 WI-4 Job A — Mpro non-covalent at 30:1 (real measured inactives)",
        "structure": "7L11 chain A", "engine": ["AutoDock Vina"] if real_vina else [],
        "real_vina": real_vina, "n_docks_timed_out": n_timeout,
        "n_actives_docked": n_act, "n_inactives_docked": len(labels) - n_act,
        "ratio": f"{(len(labels) - n_act)}:{n_act}",
        "ef1_ceiling": round(len(labels) / n_act, 1) if n_act else None,
        "benchmark": report, "bootstrap": bootstrap,
        "baseline_tables": baselines,
        "vina_vs_2d_similarity": vina_vs_2d, "holm_bonferroni": holm,
        "structure_based_enrichment_demonstrated": structure_skill,
        "plain_language": (
            "At 30:1, Vina beats 2D-similarity on BEDROC (paired-Δ CI > 0, Holm) — structure-"
            "based enrichment demonstrated." if structure_skill else
            "At 30:1, Vina still does NOT beat 2D-similarity on BEDROC — structure-based "
            "enrichment not demonstrated even at a non-degenerate ratio."),
        "note": "Separate artifact; frozen 1:1 headline + locked_benchmark.json untouched.",
        "scored_compounds": [{"name": e["name"], "dG": e["dG"],
                              "positive_control": e["positive_control"]} for e in entries],
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2))
    _md(payload)
    return payload


def _md(p: dict) -> None:
    def ci(t, m):
        r = t[m]
        return f"{r.get('median')} {r.get('ci95')}"
    b = p["benchmark"]
    md = [f"# Phase 11 WI-4 Job A — Mpro non-covalent at {p['ratio']}", "",
          f"Real Vina: {p['real_vina']}. {p['n_actives_docked']} actives / "
          f"{p['n_inactives_docked']} measured inactives. EF1% ceiling now {p['ef1_ceiling']}.", "",
          "## Baselines vs Vina (median [95% CI])", "",
          "| Method | BEDROC | logAUC | ROC-AUC | EF1% |", "|---|---|---|---|---|"]
    for label in ("random", "2d_sim", "vina"):
        t = p["baseline_tables"][label]
        md.append(f"| {label} | {ci(t,'bedroc')} | {ci(t,'logauc')} | {ci(t,'roc_auc')} | {ci(t,'ef1')} |")
    v = p["vina_vs_2d_similarity"]["bedroc"]
    md += ["", f"Paired Vina−2Dsim (BEDROC): Δ {v['median_delta']}, 95% CI {v['ci95']}, "
           f"P(Δ>0) {v['p_gt0']}.", "", f"**{p['plain_language']}**"]
    OUT.with_suffix(".md").write_text("\n".join(md) + "\n")


def main() -> None:
    p = run()
    b = p["benchmark"]
    print(f"real_vina={p['real_vina']} ratio={p['ratio']} EF1%_ceiling={p['ef1_ceiling']}")
    print(f"BEDROC={b['bedroc']} EF1%={b['ef'].get('EF1%')} ROC-AUC={b['roc_auc']}")
    print(f"structure_skill={p['structure_based_enrichment_demonstrated']}: {p['plain_language']}")
    print(OUT)


if __name__ == "__main__":
    main()
