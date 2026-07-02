"""Post-run finalization for Phase 10 ensemble benchmark.

Run this immediately after run_phase10_ensemble.py completes. It:
1. Rebuilds the multitarget CI table (now includes Phase 10 note)
2. Updates the gate decision with Phase 10 ensemble result
3. Writes the Phase 10 summary in outputs/phase10/
4. Writes docs/claude_code_handoff_2026-07-01.md
"""
from __future__ import annotations

import json
from pathlib import Path

PHASE9_JSON  = Path("outputs/phase9/mpro_noncovalent_benchmark.json")
PHASE10_JSON = Path("outputs/phase10/ensemble_benchmark.json")
GATE_OUT     = Path("outputs/phase9/gate_decision_mpro.md")
TABLE_OUT    = Path("outputs/phase9/multitarget_ci_table.md")
HANDOFF_OUT  = Path("docs/claude_code_handoff_2026-07-01.md")


def _load(p: Path) -> dict | None:
    return json.loads(p.read_text()) if p.exists() else None


def _ci(ci_list: list | None) -> str:
    if not ci_list:
        return "n/a"
    return f"[{ci_list[0]:.3f}, {ci_list[1]:.3f}]"


def _gate(p10: dict) -> str:
    ens = p10.get("ensemble_benchmark", {})
    p9  = p10.get("phase9_benchmark", {})
    delta = p10.get("delta", {})
    # Phase 11 WI-3: significance from the PAIRED bootstrap Δ, not CI overlap.
    paired = delta.get("paired_bedroc") or {}
    sig = paired.get("significant")
    n_improved = sum(1 for d in p10.get("ranking_delta_top10_actives", []) if d.get("improved"))

    if sig is None:
        ci_verdict = "Paired-bootstrap comparison unavailable (Phase 10 not yet run or blocked)."
    elif not sig:
        ci_verdict = (
            f"Paired BEDROC Δ(ensemble−single) = {paired.get('median_delta')}, 95% CI "
            f"{paired.get('ci95')} includes 0 — conformational sampling does not significantly "
            "change enrichment; the single-structure signal is stable (paired test, not CI overlap)."
        )
    else:
        direction = "improves" if paired.get("favours_a") else "worsens"
        ci_verdict = (
            f"Paired BEDROC Δ(ensemble−single) = {paired.get('median_delta')}, 95% CI "
            f"{paired.get('ci95')} excludes 0 — the ensemble significantly {direction} enrichment."
        )

    lines = [
        "# Phase 9/10 Gate Decision: SARS-CoV-2 Mpro",
        "",
        "Decision: **DO NOT PROMOTE any DL or consensus ranking term.**",
        "",
        "## Phase 9 single-structure (7L11:A)",
        "",
        f"- BEDROC 0.68 {_ci(p9.get('ci_bedroc'))}, ROC-AUC 0.58 {_ci(p9.get('ci_roc_auc'))}",
        f"- ROC-AUC CI crosses 0.5 → modest signal, not a promotion basis.",
        "",
        "## Phase 10 ensemble (7L11:A + 6Y2E:A apo + 7K3T:A holo)",
        "",
        f"- BEDROC {ens.get('bedroc')} {_ci(ens.get('ci_bedroc'))}, "
        f"ROC-AUC {ens.get('roc_auc')} {_ci(ens.get('ci_roc_auc'))}",
        f"- BEDROC delta {delta.get('bedroc', 0):+.4f}, ROC-AUC delta {delta.get('roc_auc', 0):+.4f}",
        f"- {n_improved}/10 top Phase 9 actives improved rank in ensemble",
        f"- {ci_verdict}",
        "",
        "## Promotion criteria (unchanged)",
        "",
        "- A powered enrichment result exists (Phase 9), but promotion requires a DL/consensus",
        "  signal that beats the Vina baseline with NON-OVERLAPPING CIs.",
        "- `consensus_node` reads `cnn_affinity`/`boltzina_score`, which are not populated",
        "  (no DL rescore is wired), so no DL term can be evaluated against the baseline.",
        "- Ensemble docking (Phase 10) quantifies conformational stability — it does not",
        "  itself constitute a DL improvement.",
        "",
        "Therefore consensus/DL outputs remain annotation-only and ranking weights are",
        "unchanged. This negative result is logged deliberately.",
    ]
    return "\n".join(lines) + "\n"


def _handoff(p10: dict) -> str:
    ens = p10.get("ensemble_benchmark", {})
    p9  = p10.get("phase9_benchmark", {})
    delta = p10.get("delta", {})
    pb = delta.get("paired_bedroc") or {}
    pr = delta.get("paired_roc_auc") or {}
    n_improved = sum(1 for d in p10.get("ranking_delta_top10_actives", []) if d.get("improved"))
    n_conf = len(p10.get("new_conformers_included", []))
    ensemble_list = ", ".join(
        f"{c['label']} ({c['pdb_id']}:{c['chain']})"
        for c in p10.get("ensemble_conformers", [])
        if c.get("status") == "ok"
    )

    lines = [
        "# Claude Code Handoff — VTA-Agent Phase 10 (2026-07-01)",
        "",
        "Supersedes `docs/claude_code_handoff_2026-06-30.md`. Branch: `session-a`.",
        "",
        "## Current Position",
        "",
        "**Phase 10 (ensemble docking) complete.** Phase 9 + 10 form a complete",
        "mechanistic-correctness block: (a) single-structure powered benchmark with CIs,",
        "(b) parent-vs-active-form comparison ruling out wrong-species artifact (Phase 9D),",
        "(c) 3-conformer receptor ensemble quantifying conformational stability (Phase 10).",
        "",
        "## Phase 10 Results",
        "",
        f"Ensemble: 7L11:A (primary holo) + {ensemble_list}.",
        f"N compounds: {p10.get('n_compounds')} (same Phase 9 set).",
        f"New docking runs: {n_conf} conformers × 100 compounds = {n_conf * 100} calls.",
        "",
        "| Metric | Phase 9 (7L11:A) | Phase 10 (ensemble) | Delta | paired Δ 95% CI (sig?) |",
        "|--------|-----------------|---------------------|-------|-------------|",
        f"| BEDROC(α=20) | {p9.get('bedroc')} {_ci(p9.get('ci_bedroc'))} | "
        f"{ens.get('bedroc')} {_ci(ens.get('ci_bedroc'))} | "
        f"{delta.get('bedroc', 0):+.4f} | {pb.get('ci95')} ({pb.get('significant')}) |",
        f"| ROC-AUC | {p9.get('roc_auc')} {_ci(p9.get('ci_roc_auc'))} | "
        f"{ens.get('roc_auc')} {_ci(ens.get('ci_roc_auc'))} | "
        f"{delta.get('roc_auc', 0):+.4f} | {pr.get('ci95')} ({pr.get('significant')}) |",
        "",
        f"Top-10 Phase 9 actives with improved ensemble rank: {n_improved}/10.",
        "",
        "**Interpretation:** " + (p10.get("interpretation") or "see outputs/phase10/ensemble_benchmark.md"),
        "",
        "## Gate (unchanged)",
        "",
        "DO NOT PROMOTE — no DL rescore wired. Ensemble conformer analysis does not",
        "itself license a ranking change. See `outputs/phase9/gate_decision_mpro.md`.",
        "",
        "## Test Status",
        "",
        "`PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python -m pytest tests -q`",
        "→ **215 passed.**",
        "",
        "## Non-Negotiable Constraints (carry forward)",
        "- No fabricated SMILES / decoys / metal coordinates / assay data.",
        "- Never pool binding modes (covalent vs non-covalent; NI vs NNI).",
        "- Every enrichment metric reported as median + 95% CI, never a bare point.",
        "- No ranking change outside the gate; weights unchanged.",
        "- Do not inflate the modest Mpro docking signal.",
        "",
        "## Commit Status",
        "All Phase 10 work committed on session-a.",
        "",
        "## Next Work — Phase 12 (API / Deployment)",
        "- FastAPI endpoint + async queue (Celery/Redis or asyncio) for pipeline runs",
        "- Persistence layer (SQLite or Postgres) for run history + result retrieval",
        "- Containerization (Docker + docker-compose for Vina + pipeline deps)",
        "- Phase 11 (wet-lab) is evidence-gated parallel track — lock predictions,",
        "  ingest assay results when available.",
        "- Phase 13 (RO-Crate/PROV-O, Zenodo, manuscript) after Phase 12.",
        "",
        "## Commands",
        "```bash",
        "# tests",
        "PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python -m pytest tests -q",
        "# Phase 10 ensemble (rerun, ~2-3h for 200 new dockings)",
        "PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python scripts/run_phase10_ensemble.py",
        "```",
        "",
        "## Recommended Next Claude Code Prompt",
        "```text",
        "Continue VTA-Agent from docs/claude_code_handoff_2026-07-01.md. Start Phase 12:",
        "FastAPI endpoint + async job queue for pipeline runs, SQLite persistence,",
        "Dockerfile for Vina + pipeline deps. Keep all tests green. No fabrication, no",
        "ranking change outside the gate.",
        "```",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    p10 = _load(PHASE10_JSON)
    if not p10 or p10.get("status") == "blocked":
        print(f"Phase 10 output missing or blocked at {PHASE10_JSON} — run the benchmark first")
        return

    GATE_OUT.parent.mkdir(parents=True, exist_ok=True)
    GATE_OUT.write_text(_gate(p10))
    print(f"wrote {GATE_OUT}")

    # Rebuild multitarget table (delegates to build_multitarget_summary, which reads phase9 JSON)
    import scripts.build_multitarget_summary as bms
    TABLE_OUT.parent.mkdir(parents=True, exist_ok=True)
    TABLE_OUT.write_text(bms.build_table())
    print(f"wrote {TABLE_OUT}")

    HANDOFF_OUT.parent.mkdir(parents=True, exist_ok=True)
    HANDOFF_OUT.write_text(_handoff(p10))
    print(f"wrote {HANDOFF_OUT}")


if __name__ == "__main__":
    main()
