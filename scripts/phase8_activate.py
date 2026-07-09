"""Phase 8 activation for the committed TiLV PB1 validation target."""
from __future__ import annotations

import json
from pathlib import Path

from vta.eval.decoys import registered_decoy_sets
from vta.eval.leakage import gate_recalibration_criteria, leakage_audit
from vta.eval.metrics import bootstrap_enrichment_report, enrichment_report
from vta.experiment.registry import lock_predictions

TARGET = {
    "name": "TiLV PB1",
    "structure": "8PSO chain B",
    "pocket": "NTP catalytic site",
    "eval_target_id": "TiLV_PB1_8PSO_B",
}
OUT = Path("outputs/phase8")
VALIDATION = Path("outputs/validation_controls.json")
MATCHED_BENCHMARK = OUT / "matched_benchmark_tilv_pb1.json"


def _load_validation_rows(path: Path = VALIDATION) -> list[dict]:
    if not path.exists():
        raise FileNotFoundError(f"{path} not found; run scripts/validate_controls.py first")
    rows = json.loads(path.read_text()).get("rows") or []
    if not rows:
        raise ValueError(f"{path} contains no scored rows")
    return rows


def _entries(rows: list[dict]) -> list[dict]:
    return [{
        "name": r["ligand"],
        "score": r["score"],
        "positive_control": bool(r.get("positive_control")),
        "smiles": r.get("smiles"),
        "dG": r.get("dG"),
        "species_source": r.get("species_source"),
    } for r in rows]


def _decoy_metadata(name: str = "demo") -> dict:
    meta = registered_decoy_sets().get(name)
    if not meta:
        return {"name": name, "available": False, "property_matched": False}
    return {
        "name": meta.name,
        "source": meta.source,
        "path": meta.path,
        "property_matched": meta.property_matched,
        "caveat": meta.caveat,
        "license": meta.license,
        "available": Path(meta.path).exists(),
    }


def _verdict(report: dict, decoys: dict) -> str:
    if decoys.get("property_matched") and report["n_decoys"] >= 120:
        return "publication-grade candidate"
    if decoys.get("property_matched"):
        return "defensible demonstration; underpowered for publication"
    return "demonstration only; property-matched decoys still required"


def _write_benchmark(report: dict, rows: list[dict], decoys: dict, leakage: dict) -> None:
    payload = {
        "target": TARGET,
        "benchmark": report,
        "decoy_provenance": decoys,
        "leakage_audit": leakage,
        "verdict": _verdict(report, decoys),
        "control_gate": {
            "controls_in_top": sum(1 for r in rows[:5] if r.get("positive_control")),
            "controls_total": sum(1 for r in rows if r.get("positive_control")),
            "top_n": 5,
        },
    }
    (OUT / "benchmark_tilv_pb1.json").write_text(json.dumps(payload, indent=2))
    ef = report["ef"]
    md = [
        "# Phase 8 Benchmark: TiLV PB1",
        "",
        f"Target: {TARGET['name']} ({TARGET['structure']}), {TARGET['pocket']}.",
        f"Verdict: {payload['verdict']}.",
        "",
        "## Primary Metrics",
        "",
        f"- N: {report['n']} ({report['n_actives']} actives / {report['n_decoys']} decoys)",
        f"- BEDROC(alpha=20): {report['bedroc']}",
        f"- EF1%: {ef.get('EF1%')}",
        f"- logAUC: {report['log_auc']}",
        f"- ROC-AUC: {report['roc_auc']} (secondary)",
        "",
        "## Decoy Provenance",
        "",
        f"- Source: {decoys.get('source')}",
        f"- Property matched: {decoys.get('property_matched')}",
        f"- Caveat: {decoys.get('caveat')}",
        "",
        "## Leakage Audit",
        "",
        f"- Status: {leakage['status']}",
        f"- Overlap: {leakage['overlap']}",
        "",
        "No DL scorer is promoted from this run.",
    ]
    (OUT / "benchmark_tilv_pb1.md").write_text("\n".join(md) + "\n")


def _write_predictions(rows: list[dict]) -> dict:
    registry = lock_predictions(rows, "phase8_tilv_pb1", top_n=min(10, len(rows)))
    registry["target"] = TARGET
    (OUT / "locked_predictions_tilv_pb1.json").write_text(json.dumps(registry, indent=2))
    return registry


def _write_experimental_plan(registry: dict) -> None:
    top = registry["predictions"][0]["ligand"] if registry["predictions"] else "top lead"
    md = [
        "# Phase 8 Prospective Experimental Plan: TiLV PB1",
        "",
        f"Locked registry: `{registry['run_id']}` at {registry['locked_at']}.",
        "",
        "## Minimal Assay Package",
        "",
        "- Primary: biochemical TiLV PB1/RdRp inhibition IC50 if assay is available.",
        "- Alternative: cell antiviral EC50 plus CC50 to compute selectivity index.",
        "- Orthogonal binding confirmation: thermal shift, SPR, or ITC when feasible.",
        "",
        "## Success Criterion",
        "",
        f"- At least one locked top prediction, starting with {top}, shows reproducible",
        "  activity in the selected assay with a measurable dose-response.",
        "- For cell assays, require CC50/EC50 selectivity index above 10 for a credible hit.",
        "",
        "## Ownership",
        "",
        "- Assay owner/collaborator: TBD.",
        "- No retrospective edits to the locked prediction set are allowed.",
    ]
    (OUT / "experimental_plan_tilv_pb1.md").write_text("\n".join(md) + "\n")


def _write_gate_decision(report: dict, decoys: dict, leakage: dict) -> None:
    criteria = gate_recalibration_criteria()
    promote = bool(decoys.get("property_matched") and leakage["promotable"] and report["n_decoys"] >= 120)
    md = [
        "# Phase 8 Gate Decision: DL Promotion",
        "",
        f"Decision: {'PROMOTE candidate term' if promote else 'DO NOT PROMOTE any DL term'}.",
        "",
        "Reason:",
        "",
        f"- Property-matched decoys: {decoys.get('property_matched')}",
        f"- Leakage audit promotable: {leakage['promotable']}",
        f"- Decoy count: {report['n_decoys']}",
        f"- Required gate items: {', '.join(criteria['required'])}",
        "",
        "Ranking remains unchanged. Consensus/DL outputs stay annotation-only.",
    ]
    (OUT / "gate_decision_tilv_pb1.md").write_text("\n".join(md) + "\n")


def _write_matched_summary_if_available() -> bool:
    if not MATCHED_BENCHMARK.exists():
        return False
    payload = json.loads(MATCHED_BENCHMARK.read_text())
    bench = payload["benchmark"]
    ci = payload.get("bootstrap", {}).get("bootstrap") or {}
    leakage = payload["leakage_audit"]
    controls_top5 = sum(1 for r in payload["rows"][:5] if r.get("positive_control"))
    controls_total = sum(1 for r in payload["rows"] if r.get("positive_control"))
    md = [
        "# Phase 8 Matched Benchmark Summary",
        "",
        f"Target: {payload['target']}.",
        f"Verdict: {payload['verdict']}.",
        "",
        "## Primary Metrics",
        "",
        f"- N: {bench['n']} ({bench['n_actives']} actives / {bench['n_decoys']} decoys)",
        f"- BEDROC(alpha=20): {bench['bedroc']}",
        f"- EF1%: {bench['ef'].get('EF1%')}",
        f"- logAUC: {bench['log_auc']}",
        f"- ROC-AUC: {bench['roc_auc']} (secondary)",
        f"- Controls in top 5: {controls_top5}/{controls_total}",
        "",
        "## Leakage Audit",
        "",
        f"- Status: {leakage['status']}",
        f"- Overlap: {leakage['overlap']}",
        "",
        "This matched benchmark is the current Phase 8 headline artifact. It is still",
        "underpowered and does not justify a ranking change.",
    ]
    (OUT / "benchmark_matched_tilv_pb1.md").write_text("\n".join(md) + "\n")
    if ci:
        _write_ci_summary(ci)
    _write_matched_gate_decision(payload)
    return True


def _write_ci_summary(ci: dict) -> None:
    lines = ["# Phase 9 Bootstrap CI Summary", ""]
    for metric in ("bedroc", "EF1%", "log_auc", "roc_auc"):
        row = ci.get(metric) or {}
        lines.append(
            f"- {metric}: median {row.get('median')}, 95% CI {row.get('ci95')}")
    lines.append("")
    lines.append("Intervals are bootstrap estimates over compounds and remain unstable")
    lines.append("until the active set is expanded beyond the current four positives.")
    (OUT / "benchmark_ci_tilv_pb1.md").write_text("\n".join(lines) + "\n")


def _write_matched_gate_decision(payload: dict) -> None:
    bench = payload["benchmark"]
    leakage = payload["leakage_audit"]
    promote = bool(leakage["promotable"] and bench["n_decoys"] >= 120)
    md = [
        "# Phase 8 Gate Decision: Matched Benchmark",
        "",
        f"Decision: {'PROMOTE candidate term' if promote else 'DO NOT PROMOTE any DL or consensus term'}.",
        "",
        "Evidence:",
        "",
        f"- Target: {payload['target']}.",
        "- Benchmark file: `outputs/phase8/matched_benchmark_tilv_pb1.json`",
        f"- N = {bench['n']}: {bench['n_actives']} actives / {bench['n_decoys']} decoys.",
        f"- BEDROC(alpha=20) = {bench['bedroc']}",
        f"- EF1% = {bench['ef'].get('EF1%')}",
        f"- logAUC = {bench['log_auc']}",
        f"- ROC-AUC = {bench['roc_auc']}",
        "",
        "Reason:",
        "",
        "- The benchmark is underpowered for ranking-term promotion.",
        "- Decoys are presumed inactive, not experimentally verified inactive.",
        "- GNINA/Boltz-2 training-target exclusion data were not supplied.",
        "- No leakage-controlled DL-vs-baseline enrichment improvement has been shown.",
        "",
        "Ranking remains unchanged.",
    ]
    (OUT / "gate_decision_matched_tilv_pb1.md").write_text("\n".join(md) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = sorted(_load_validation_rows(), key=lambda r: r["score"], reverse=True)
    report = enrichment_report(_entries(rows))
    ci_report = bootstrap_enrichment_report(_entries(rows), n_resamples=1000, seed=42)
    decoys = _decoy_metadata("demo")
    leakage = leakage_audit([], [TARGET["eval_target_id"]])
    _write_benchmark(report, rows, decoys, leakage)
    (OUT / "benchmark_tilv_pb1_ci.json").write_text(json.dumps(ci_report, indent=2))
    registry = _write_predictions(rows)
    _write_experimental_plan(registry)
    _write_gate_decision(report, decoys, leakage)
    has_matched = _write_matched_summary_if_available()
    print(f"Phase 8 target: {TARGET['name']} / {TARGET['structure']}")
    print(f"Benchmark verdict: {_verdict(report, decoys)}")
    print(f"BEDROC(alpha=20)={report['bedroc']} EF1%={report['ef'].get('EF1%')} logAUC={report['log_auc']}")
    print(f"Locked predictions: {len(registry['predictions'])} -> {OUT / 'locked_predictions_tilv_pb1.json'}")
    if has_matched:
        print(f"Matched benchmark summary -> {OUT / 'benchmark_matched_tilv_pb1.md'}")
    print(f"Wrote artifacts to {OUT}")


if __name__ == "__main__":
    main()
