"""Powered SARS-CoV-2 Mpro enrichment benchmark (Phase 9C, Track B).

Real AutoDock Vina over **non-covalent** Moonshot actives plus **experimentally measured**
Moonshot inactives, docked into the 7L11 catalytic active site. This is the powered,
CI-informative counterpart to the underpowered TiLV PB1 demonstration and the
un-benchmarkable HCV nucleotide set.

Why Moonshot-only for the headline: actives and inactives come from the *same* fluorescence
assay, so the active/inactive label is not confounded by cross-assay chemistry. ChEMBL
(CHEMBL4523582) informed the dataset and is reported as a cross-check; it is only used to
top up actives if Moonshot falls short (it does not here).

Ranking signal = Vina affinity (score = −ΔG): the cleanest docking-enrichment number, not
confounded by the placeholder conservation term. Metrics are reported as median + 95% CI
(bootstrap over compounds) plus a multi-draw control-spread across inactive subsamples.
"""
from __future__ import annotations

import argparse
import json
import random
import statistics
from pathlib import Path

import vta.nodes.docking as docking
from vta.eval.leakage import leakage_audit
from vta.eval.metrics import bootstrap_enrichment_report, enrichment_report
from vta.nodes.docking import docking_node
from vta.nodes.pockets import pockets_node
from vta.nodes.rank import rank_node
from vta.nodes.structure import structure_node
from vta.state import new_state

STRAT = Path("vta/data/mpro/mpro_stratified.json")
OUT = Path("outputs/phase9/mpro_noncovalent_benchmark.json")
EVAL_TARGET_ID = "MPRO_7L11_A"


def _sanitize(text: str) -> str:
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in (text or ""))[:48]


def select_library(strat: dict, n_actives: int, n_inactives: int, seed: int,
                   prefer_source: str = "moonshot") -> tuple[list[dict], list[dict]]:
    """Seeded random sample of non-covalent actives + measured inactives (no cherry-pick)."""
    actives = [r for r in strat["non_covalent_actives"] if r["source"] == prefer_source]
    inactives = [r for r in strat["inactives"] if r["source"] == prefer_source]
    if len(actives) < n_actives:  # top up actives from ChEMBL only if Moonshot is short
        extra = [r for r in strat["non_covalent_actives"] if r["source"] != prefer_source]
        random.Random(seed).shuffle(extra)
        actives = actives + extra
    rng = random.Random(seed)
    rng.shuffle(actives)
    rng.shuffle(inactives)
    return actives[:n_actives], inactives[:n_inactives]


def to_ligand(rec: dict, is_active: bool, idx: int) -> dict:
    sid = _sanitize(rec.get("id") or "") or f"MPRO_{idx}"
    return {
        "name": rec.get("id") or sid,
        "chembl_id": sid,
        "dock_cache_id": sid,
        "smiles": rec["smiles"],
        "max_phase": None,
        "positive_control": is_active,
        "source": rec["source"],
    }


def _entries(rows: list[dict]) -> list[dict]:
    """Enrichment entries scored by Vina affinity (−ΔG: higher = more active-like)."""
    out = []
    for r in rows:
        dG = r.get("dG")
        if dG is None:
            continue
        out.append({
            "name": r["ligand"],
            "score": round(-float(dG), 4),
            "positive_control": bool(r.get("positive_control")),
            "smiles": r.get("smiles"),
            "dG": dG,
            "composite_rank_score": r.get("score"),
        })
    return out


def control_draw_spread(entries: list[dict], *, draws: int = 6, seed: int = 7,
                        keep_fraction: float = 0.8, alpha: float = 20.0) -> dict:
    """Spread of metrics across random inactive subsamples (robustness to control choice).

    Keeps all actives; for each draw, samples `keep_fraction` of the inactives without
    replacement and recomputes the metrics. Reports min/median/max across draws.
    """
    actives = [e for e in entries if e.get("positive_control")]
    inactives = [e for e in entries if not e.get("positive_control")]
    if not actives or not inactives:
        return {"draws": 0, "note": "need both classes"}
    rng = random.Random(seed)
    keep = max(2, int(round(keep_fraction * len(inactives))))
    collected: dict[str, list[float]] = {"bedroc": [], "EF1%": [], "log_auc": [], "roc_auc": []}
    for _ in range(draws):
        sample_inactives = rng.sample(inactives, min(keep, len(inactives)))
        rep = enrichment_report(actives + sample_inactives, fractions=(0.01,), alpha=alpha)
        collected["bedroc"].append(rep["bedroc"])
        collected["EF1%"].append(rep["ef"]["EF1%"])
        collected["log_auc"].append(rep["log_auc"])
        collected["roc_auc"].append(rep["roc_auc"])
    spread = {
        name: {
            "min": round(min(vals), 4),
            "median": round(statistics.median(vals), 4),
            "max": round(max(vals), 4),
        }
        for name, vals in collected.items()
    }
    return {"draws": draws, "keep_fraction": keep_fraction, "n_inactives_per_draw": keep,
            "seed": seed, "spread": spread}


def run(n_actives: int = 50, n_inactives: int = 50, seed: int = 17,
        n_resamples: int = 2000) -> dict:
    strat = json.loads(STRAT.read_text())
    sel_actives, sel_inactives = select_library(strat, n_actives, n_inactives, seed)
    library = ([to_ligand(r, True, i) for i, r in enumerate(sel_actives)]
               + [to_ligand(r, False, i + 10000) for i, r in enumerate(sel_inactives)])
    docking.load_ligands = lambda: library

    st = new_state("mpro_noncovalent_benchmark", "mpro_noncovalent_benchmark")
    st["extracted_proteins"] = {"MPRO": {"sequence": "M" * 306, "length_aa": 306, "plddt": None}}
    st = structure_node(st)
    st = pockets_node(st)
    st = docking_node(st)
    st = rank_node(st)

    rows = sorted(st["docking_results"], key=lambda r: (-(r["dG"]) if r.get("dG") is not None else -1e9),
                  reverse=True)
    entries = _entries(rows)
    engines = sorted({(r.get("dG_provenance") or {}).get("engine") for r in rows})
    real_vina = engines == ["AutoDock Vina"]

    # Honesty guard: if no real Vina pose was produced (e.g. the known Meeko 0.7.1 Mpro
    # receptor-prep failure), do NOT emit a fake/mock enrichment number. Write a blocked
    # status artifact instead, recording what is ready and the exact blocker.
    if not entries or not real_vina:
        return _write_blocked(st, sel_actives, sel_inactives, n_actives, n_inactives, seed)
    report = enrichment_report(entries)
    bootstrap = bootstrap_enrichment_report(entries, n_resamples=n_resamples, seed=42)
    spread = control_draw_spread(entries)

    n_actives_docked = sum(1 for e in entries if e["positive_control"])
    n_inactives_docked = sum(1 for e in entries if not e["positive_control"])
    ci_bedroc = (bootstrap.get("bootstrap", {}).get("bedroc", {}) or {}).get("ci95")
    informative = bool(ci_bedroc and (ci_bedroc[1] - ci_bedroc[0]) < 0.99)

    payload = {
        "target": "SARS-CoV-2 Mpro (3CLpro), non-covalent",
        "structure": "7L11 chain A",
        "active_site": strat.get("headline_benchmark_inputs"),
        "binding_mode": "non_covalent actives only (covalent held out)",
        "control_set": "experimentally measured Moonshot inactives (real, not presumed decoys)",
        "selection": {"requested_actives": n_actives, "requested_inactives": n_inactives,
                      "seed": seed, "prefer_source": "moonshot"},
        "n_actives_selected": len(sel_actives),
        "n_inactives_selected": len(sel_inactives),
        "n_actives_docked": n_actives_docked,
        "n_inactives_docked": n_inactives_docked,
        "engine": engines,
        "real_vina": real_vina,
        "scoring": "enrichment ranked by Vina affinity (score = -dG)",
        "benchmark": report,
        "bootstrap": bootstrap,
        "control_draw_spread": spread,
        "ci_informative": informative,
        "leakage_audit": leakage_audit([], [EVAL_TARGET_ID]),
        "versions": st.get("versions"),
        "rows": rows,
        "verdict": (
            "powered, publication-grade candidate (real measured inactives, informative CIs)"
            if (real_vina and informative and n_actives_docked >= 30)
            else "demonstration-grade (check engine/CI/active count)"
        ),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2))
    _write_md(payload)
    return payload


def _write_blocked(st: dict, sel_actives: list[dict], sel_inactives: list[dict],
                   n_actives: int, n_inactives: int, seed: int) -> dict:
    """Honest 'blocked' artifact: dataset is ready, real docking is gated on receptor prep."""
    prep_notes = [a for a in (st.get("audit_trail") or []) if "Vina[MPRO]" in a or "receptor" in a]
    payload = {
        "target": "SARS-CoV-2 Mpro (3CLpro), non-covalent",
        "structure": "7L11 chain A",
        "status": "blocked_receptor_prep",
        "blocker": (
            "AutoDock Vina is installed, but Meeko 0.7.1 receptor preparation fails "
            "deterministically on the Mpro chain (mk_prepare_receptor raises "
            "'update_H_positions: Updated 1 H positions but deleted N' for every Mpro PDB "
            "tried: 6Y2E, 6M03, 7K3T, 7TLL, 7L10, 6W63, 7RFW, 7VH8, 7L11, 5R8T, 7BB2). "
            "The same code path preps the TiLV 8PSO receptor successfully, so the engine and "
            "pipeline are intact; the failure is Meeko-vs-Mpro specific. No fallback receptor "
            "prep tool (OpenBabel / reduce / ADFR) is installed."
        ),
        "what_is_ready": {
            "dataset": "vta/data/mpro/mpro_dataset.json",
            "stratified": "vta/data/mpro/mpro_stratified.json",
            "active_site_wiring": "EXPERIMENTAL_PDB/EXPERIMENTAL_ACTIVE_SITE += MPRO (7L11:A)",
            "selected_actives": len(sel_actives),
            "selected_inactives": len(sel_inactives),
            "binding_mode": "non_covalent actives only (covalent held out)",
            "control_set": "experimentally measured Moonshot inactives",
            "selection_seed": seed,
        },
        "next_step": (
            "Provide a working receptor-prep path (OpenBabel fallback seam, or a Meeko "
            "version without this bug), then re-run scripts/run_mpro_benchmark.py to produce "
            "the powered, CI-bearing enrichment number. No mock/fabricated number is emitted."
        ),
        "receptor_prep_audit": prep_notes,
        "verdict": "blocked — powered benchmark not yet runnable; no number claimed",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2))
    md = [
        "# Phase 9C Benchmark: SARS-CoV-2 Mpro (non-covalent) — BLOCKED",
        "",
        f"Status: **{payload['status']}**.",
        f"Verdict: {payload['verdict']}.",
        "",
        "## Blocker",
        "",
        payload["blocker"],
        "",
        "## What is ready (no docking required)",
        "",
        f"- Dataset: `{payload['what_is_ready']['dataset']}`",
        f"- Stratified (binding mode): `{payload['what_is_ready']['stratified']}`",
        f"- Selected non-covalent actives: {payload['what_is_ready']['selected_actives']}",
        f"- Selected measured inactives: {payload['what_is_ready']['selected_inactives']}",
        f"- Active-site wiring: {payload['what_is_ready']['active_site_wiring']}",
        "",
        "## Next step",
        "",
        payload["next_step"],
        "",
        "No mock or fabricated enrichment number is reported for Mpro.",
    ]
    OUT.with_suffix(".md").write_text("\n".join(md) + "\n")
    return payload


def _write_md(payload: dict) -> None:
    bench = payload["benchmark"]
    boot = payload.get("bootstrap", {}).get("bootstrap", {}) or {}

    def ci(metric: str) -> str:
        row = boot.get(metric) or {}
        return f"median {row.get('median')}, 95% CI {row.get('ci95')}" if row else "n/a"

    md = [
        "# Phase 9C Benchmark: SARS-CoV-2 Mpro (non-covalent)",
        "",
        f"Target: {payload['target']}. Structure: {payload['structure']}.",
        f"Binding mode: {payload['binding_mode']}.",
        f"Control set: {payload['control_set']}.",
        f"Engine: {', '.join(e for e in payload['engine'] if e)} "
        f"(real Vina: {payload['real_vina']}).",
        f"Verdict: {payload['verdict']}.",
        "",
        "## Primary metrics (point + bootstrap CI)",
        "",
        f"- N: {bench['n']} ({bench['n_actives']} actives / {bench['n_decoys']} inactives)",
        f"- BEDROC(alpha=20): {bench['bedroc']} — {ci('bedroc')}",
        f"- EF1%: {bench['ef'].get('EF1%')} — {ci('EF1%')}",
        f"- logAUC: {bench['log_auc']} — {ci('log_auc')}",
        f"- ROC-AUC: {bench['roc_auc']} — {ci('roc_auc')} (secondary)",
        f"- CIs informative (not [0,1]): {payload['ci_informative']}",
        "",
        "## Control-draw spread (robustness to inactive subsample)",
        "",
    ]
    spread = (payload.get("control_draw_spread") or {}).get("spread") or {}
    for metric, vals in spread.items():
        md.append(f"- {metric}: median {vals['median']} (min {vals['min']}, max {vals['max']})")
    md += [
        "",
        "## Leakage audit",
        "",
        f"- Status: {payload['leakage_audit']['status']}",
        f"- Overlap: {payload['leakage_audit']['overlap']}",
        "",
        "Ranking weights unchanged. Consensus/DL outputs remain annotation-only.",
    ]
    OUT.with_suffix(".md").write_text("\n".join(md) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-actives", type=int, default=50)
    parser.add_argument("--n-inactives", type=int, default=50)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--resamples", type=int, default=2000)
    args = parser.parse_args()
    out = run(args.n_actives, args.n_inactives, args.seed, args.resamples)
    if out.get("status") == "blocked_receptor_prep":
        print(f"STATUS: {out['status']} — {out['verdict']}")
        print(f"ready: {out['what_is_ready']['selected_actives']} actives / "
              f"{out['what_is_ready']['selected_inactives']} inactives selected")
        print(OUT)
        return
    b = out["benchmark"]
    print(f"engine={out['engine']} real_vina={out['real_vina']}")
    print(f"N={b['n']} actives={b['n_actives']} inactives={b['n_decoys']}")
    print(f"BEDROC={b['bedroc']} EF1%={b['ef'].get('EF1%')} logAUC={b['log_auc']} ROC-AUC={b['roc_auc']}")
    print(f"ci_informative={out['ci_informative']} verdict={out['verdict']}")
    print(OUT)


if __name__ == "__main__":
    main()
