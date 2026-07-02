"""Phase 11 WI-3b — re-test Phase-10 ensemble vs single-structure with the PAIRED bootstrap.

Phase 10 decided "ensemble ≈ single (CIs overlap)". That is the invalid CI-overlap test. Both
rankings score the SAME 100 compounds (ensemble ΔG vs single-structure/primary ΔG per
compound, committed in `ensemble_rows`), so the correct test is a paired bootstrap of the
per-compound difference. This reads the committed rows — no re-docking — and writes the
corrected verdict that supersedes the CI-overlap statement.
"""
from __future__ import annotations

import json
from pathlib import Path

from vta.eval.metrics import METRIC_NAMES
from vta.eval.significance import holm_bonferroni, paired_bootstrap_delta

SRC = Path("outputs/phase10/ensemble_benchmark.json")
OUT = Path("outputs/phase11/ensemble_paired.json")
N = 10000
PRIMARY = "bedroc"


def run() -> dict:
    rows = json.loads(SRC.read_text())["ensemble_rows"]
    labels = [1 if r.get("positive_control") else 0 for r in rows]
    # higher score = more active-like (−ΔG). primary_dG = single-structure 7L11; ensemble_dG = best-of-ensemble.
    ensemble = [-float(r["ensemble_dG"]) for r in rows]
    single = [-float(r.get("primary_dG", r["dG"])) for r in rows]

    fam = {m: paired_bootstrap_delta(ensemble, single, labels, m, n=N, seed=41) for m in METRIC_NAMES}
    holm = holm_bonferroni({m: (1.0 - fam[m]["p_gt0"]) if fam[m]["p_gt0"] is not None else None
                            for m in METRIC_NAMES})
    prim = fam[PRIMARY]
    significant = bool(prim["significant"] and holm[PRIMARY]["reject"])
    verdict = (
        f"Ensemble differs from single-structure on {PRIMARY} (paired Δ={prim['median_delta']}, "
        f"95% CI {prim['ci95']}, {'favours ensemble' if prim['favours_a'] else 'favours single'})."
        if significant else
        f"No significant ensemble effect on {PRIMARY}: paired Δ={prim['median_delta']}, 95% CI "
        f"{prim['ci95']} includes 0 (P(Δ>0)={prim['p_gt0']}). Single-structure signal is stable "
        f"under conformational sampling — established by a PAIRED test, not CI overlap.")
    payload = {
        "phase": "11 WI-3b — ensemble vs single, paired bootstrap",
        "supersedes": "the CI-overlap significance statement in outputs/phase10/",
        "n_compounds": len(labels), "n_actives": sum(labels),
        "primary_metric": PRIMARY, "per_metric_paired_delta": fam,
        "holm_bonferroni": holm, "significant_primary": significant, "verdict": verdict,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2))
    md = ["# Phase 11 WI-3b — Ensemble vs single-structure (paired bootstrap)", "",
          payload["verdict"], "",
          "| Metric | median Δ(ens−single) | 95% CI | P(Δ>0) | Holm reject |",
          "|---|---|---|---|---|"]
    for m in METRIC_NAMES:
        f = fam[m]
        md.append(f"| {m} | {f['median_delta']} | {f['ci95']} | {f['p_gt0']} | "
                  f"{holm[m]['reject']} |")
    md += ["", "This replaces the invalid CI-overlap comparison (Schenker & Gentleman 2001; "
           "Cumming 2009)."]
    OUT.with_suffix(".md").write_text("\n".join(md) + "\n")
    return payload


def main() -> None:
    p = run()
    print(p["verdict"])
    print(OUT)


if __name__ == "__main__":
    main()
