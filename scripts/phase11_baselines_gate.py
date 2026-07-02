"""Phase 11 WI-1/WI-3/WI-6 — trivial baselines, paired-bootstrap significance, LE gate.

Runs entirely on the committed docked benchmarks (no re-docking). For each benchmark it scores
compounds three ways where relevant — Vina (−ΔG), the composite rank, ΔG-only, LE-only — plus
the two trivial baselines (2D-similarity, random), all on the IDENTICAL actives/decoys split,
and compares them with the paired bootstrap in `vta.eval.significance` (never CI overlap).

Outputs `outputs/phase11/baselines_and_gate.{json,md}`:
  * WI-1: per-benchmark {random, 2D-sim, Vina} metric table with 95% CIs.
  * WI-3: paired Vina-vs-2D-sim and Vina-vs-random (Holm–Bonferroni across the metric family);
    plain-language "structure-based enrichment NOT demonstrated" when Vina fails to beat 2D-sim.
  * WI-6: paired composite-vs-ΔG-only and LE-only-vs-ΔG-only; a gate decision on whether the
    LE-weighted composite earns its place over plain ΔG ranking.
"""
from __future__ import annotations

import json
from pathlib import Path

from vta.eval.baselines import baseline_2d_similarity, bootstrap_metric_ci, random_metric_distribution
from vta.eval.metrics import METRIC_NAMES
from vta.eval.significance import holm_bonferroni, paired_bootstrap_delta
from vta.nodes.rank import WEIGHTS as RANK_WEIGHTS

OUT = Path("outputs/phase11")
N = 10000
PRIMARY = "bedroc"
BENCHMARKS = {
    "Mpro_noncovalent": {"path": "outputs/phase9/mpro_noncovalent_benchmark.json", "powered": True},
    "TiLV_PB1": {"path": "outputs/phase8/matched_benchmark_tilv_pb1.json", "powered": False},
}


def _rows(path: str) -> list[dict] | None:
    p = Path(path)
    if not p.exists():
        return None
    return json.loads(p.read_text()).get("rows")


def _methods(rows: list[dict]) -> tuple[dict[str, list[float]], list[int], list[str]]:
    """Aligned per-compound scores for each method (higher = more active-like)."""
    labels = [1 if r.get("positive_control") else 0 for r in rows]
    smiles = [r.get("smiles") or r.get("dock_smiles") or "" for r in rows]
    methods = {
        "vina": [-float(r["dG"]) for r in rows],                    # −ΔG
        "composite": [float(r.get("score", 0.0)) for r in rows],    # rank_node composite
        "dG_only": [-float(r["dG"]) for r in rows],
        "le_only": [-float(r["le"]) for r in rows if r.get("le") is not None]
        if all(r.get("le") is not None for r in rows) else None,
    }
    methods["2d_sim"] = baseline_2d_similarity(smiles, labels)
    return methods, labels, smiles


def _paired(a, b, labels, seed):
    fam = {m: paired_bootstrap_delta(a, b, labels, m, n=N, seed=seed) for m in METRIC_NAMES}
    holm = holm_bonferroni({m: (1.0 - fam[m]["p_gt0"]) if fam[m]["p_gt0"] is not None else None
                            for m in METRIC_NAMES})
    return {"per_metric": fam, "holm_bonferroni": holm,
            "primary_favours_a": bool(fam[PRIMARY]["favours_a"] and holm[PRIMARY]["reject"])}


def analyse(name: str, rows: list[dict], powered: bool) -> dict:
    methods, labels, _ = _methods(rows)
    n_act = sum(labels)
    # WI-1 metric tables (bootstrap CI) for the reportable methods.
    tables = {}
    for label in ("random", "2d_sim", "vina"):
        if label == "random":
            tables[label] = {m: random_metric_distribution(labels, m, n=N, seed=11) for m in METRIC_NAMES}
        else:
            tables[label] = {m: bootstrap_metric_ci(methods[label], labels, m, n=N, seed=11)
                             for m in METRIC_NAMES}
    # WI-3 paired: Vina vs 2D-sim (the memorization test) and Vina vs random.
    vina_vs_2d = _paired(methods["vina"], methods["2d_sim"], labels, seed=21)
    vina_vs_rand_note = ("random has no fixed per-compound score; see the random metric table "
                         "and Vina's CI (both scored on the same split).")
    # WI-6 gate: composite vs ΔG-only; LE-only vs ΔG-only.
    wi6 = None
    if methods["le_only"] is not None:
        composite_vs_dg = _paired(methods["composite"], methods["dG_only"], labels, seed=31)
        le_vs_dg = _paired(methods["le_only"], methods["dG_only"], labels, seed=32)
        composite_earns_le = composite_vs_dg["primary_favours_a"]
        wi6 = {
            "composite_vs_dG_only": composite_vs_dg,
            "le_only_vs_dG_only": le_vs_dg,
            "gate_decision": ("KEEP composite (LE earns its weight)" if composite_earns_le
                              else "DEMOTE LE to secondary annotation; make ΔG the primary ranking term"),
            # Truthful: reflects the live rank_node weights, not a hardcoded flag.
            "ranking_change_applied": bool(not composite_earns_le and RANK_WEIGHTS.get("le", 0) == 0
                                           and RANK_WEIGHTS.get("dG", 0) > 0),
            "live_rank_weights": dict(RANK_WEIGHTS),
            "note": ("Gated per WI-6: composite is kept only if it beats ΔG-only on BEDROC with a "
                     "paired-Δ 95% CI strictly > 0 (after Holm). Otherwise LE is demoted."),
        }
    structure_skill = vina_vs_2d["primary_favours_a"]
    return {
        "benchmark": name,
        "powered": powered,
        "n_compounds": len(labels),
        "n_actives": n_act,
        "n_inactives": len(labels) - n_act,
        "metric_family": list(METRIC_NAMES),
        "primary_metric": PRIMARY,
        "wi1_baseline_tables": tables,
        "wi3_vina_vs_2d_similarity": vina_vs_2d,
        "wi3_vina_vs_random_note": vina_vs_rand_note,
        "structure_based_enrichment_demonstrated": bool(structure_skill),
        "wi3_plain_language": (
            "structure-based enrichment demonstrated: Vina beats 2D-similarity on BEDROC "
            "(paired-Δ CI > 0 after Holm)." if structure_skill else
            "structure-based enrichment NOT demonstrated on this target: Vina does not beat "
            "2D-similarity memorization on BEDROC with a positive paired-difference CI."),
        "wi6_le_gate": wi6,
    }


def _fmt(ci_obj: dict) -> str:
    return f"{ci_obj.get('median')} {ci_obj.get('ci95')}"


def _md(results: list[dict]) -> str:
    L = ["# Phase 11 — Trivial baselines, paired significance, LE gate", "",
         f"Every number is a bootstrap median + 95% CI (n={N}). Significance is by paired "
         "bootstrap of the difference on the SAME split — never CI overlap.", ""]
    for r in results:
        L += [f"## {r['benchmark']} ({'powered' if r['powered'] else 'UNDERPOWERED'}; "
              f"{r['n_actives']} actives / {r['n_inactives']} inactives)", "",
              "### WI-1 baselines vs Vina (median [95% CI])", "",
              "| Method | BEDROC | logAUC | ROC-AUC | EF1% |", "|---|---|---|---|---|"]
        for label in ("random", "2d_sim", "vina"):
            t = r["wi1_baseline_tables"][label]
            L.append(f"| {label} | {_fmt(t['bedroc'])} | {_fmt(t['logauc'])} | "
                     f"{_fmt(t['roc_auc'])} | {_fmt(t['ef1'])} |")
        v2 = r["wi3_vina_vs_2d_similarity"]["per_metric"][PRIMARY]
        L += ["", "### WI-3 paired test — Vina vs 2D-similarity (primary: BEDROC)", "",
              f"- median Δ(Vina−2Dsim) = {v2['median_delta']}, 95% CI {v2['ci95']}, "
              f"P(Δ>0) = {v2['p_gt0']}",
              f"- **{r['wi3_plain_language']}**", ""]
        if r["wi6_le_gate"]:
            g = r["wi6_le_gate"]
            cd = g["composite_vs_dG_only"]["per_metric"][PRIMARY]
            L += ["### WI-6 LE gate (primary: BEDROC)", "",
                  f"- composite vs ΔG-only: median Δ = {cd['median_delta']}, CI {cd['ci95']}, "
                  f"P(Δ>0) = {cd['p_gt0']}",
                  f"- **Gate: {g['gate_decision']}** (ranking change applied: "
                  f"{g['ranking_change_applied']})", ""]
    return "\n".join(L) + "\n"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    results = []
    for name, meta in BENCHMARKS.items():
        rows = _rows(meta["path"])
        if not rows:
            print(f"skip {name}: no rows at {meta['path']}")
            continue
        if not all(r.get("smiles") or r.get("dock_smiles") for r in rows):
            print(f"skip {name}: rows lack SMILES (2D-sim baseline needs them)")
            continue
        results.append(analyse(name, rows, meta["powered"]))
    (OUT / "baselines_and_gate.json").write_text(json.dumps(results, indent=2))
    (OUT / "baselines_and_gate.md").write_text(_md(results))
    for r in results:
        print(f"{r['benchmark']}: structure_skill={r['structure_based_enrichment_demonstrated']} "
              f"| LE gate={r['wi6_le_gate']['gate_decision'] if r['wi6_le_gate'] else 'n/a'}")
    print(OUT / "baselines_and_gate.md")


if __name__ == "__main__":
    main()
