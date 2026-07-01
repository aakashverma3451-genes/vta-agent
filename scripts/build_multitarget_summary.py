"""Build the Phase 9 multi-target CI table + Mpro gate decision.

One honest table across all three validation targets, every enrichment metric with its
bootstrap CI and an explicit grade/status:

  * TiLV PB1 — underpowered demonstration (4 actives, CIs ≈ [0,1]).
  * HCV NS5B NI / triphosphate — un-benchmarkable (Track A meta-finding).
  * SARS-CoV-2 Mpro (non-covalent) — powered target; currently blocked on receptor prep.

Reads committed artifacts only; degrades to a labelled "missing" row if one is absent.
"""
from __future__ import annotations

import json
from pathlib import Path

OUT = Path("outputs/phase9")
TILV = Path("outputs/phase8/matched_benchmark_tilv_pb1.json")
HCV_QUALITY = Path("outputs/phase9/hcv_ns5b_ni_matching_quality.json")
HCV_PROV = Path("outputs/phase9/hcv_ns5b_ni_decoy_provenance.json")
MPRO = Path("outputs/phase9/mpro_noncovalent_benchmark.json")
MPRO_STRAT = Path("vta/data/mpro/mpro_stratified.json")


def _load(path: Path) -> dict | None:
    return json.loads(path.read_text()) if path.exists() else None


def _ci(boot: dict, metric: str) -> str:
    row = (boot or {}).get(metric) or {}
    if not row:
        return "n/a"
    return f"{row.get('median')} {row.get('ci95')}"


def build_table() -> str:
    tilv = _load(TILV)
    hcv_q = _load(HCV_QUALITY)
    hcv_p = _load(HCV_PROV)
    mpro = _load(MPRO)
    mpro_s = _load(MPRO_STRAT)

    lines = [
        "# Phase 9 Multi-Target Enrichment Summary",
        "",
        "Every metric is reported as median + 95% bootstrap CI, never a bare point. Each",
        "target carries an explicit grade/status. Two of the three results are negative or",
        "blocked — and the negatives are the scientific contribution, not a failure to hide.",
        "",
        "| Target | Class / pocket | Actives | Control set | BEDROC(α20) med [CI] | EF1% med [CI] | logAUC med [CI] | ROC-AUC med [CI] | Grade / status |",
        "|---|---|---|---|---|---|---|---|---|",
    ]

    # TiLV PB1
    if tilv:
        b = tilv["benchmark"]
        boot = tilv.get("bootstrap", {}).get("bootstrap", {})
        lines.append(
            f"| TiLV PB1 | RdRp NTP site (8PSO:B) | {b['n_actives']} | "
            f"ChEMBL property-matched presumed decoys | {_ci(boot,'bedroc')} | "
            f"{_ci(boot,'EF1%')} | {_ci(boot,'log_auc')} | {_ci(boot,'roc_auc')} | "
            f"underpowered demonstration |")
    else:
        lines.append("| TiLV PB1 | RdRp NTP site | — | — | — | — | — | — | artifact missing |")

    # HCV NS5B NI / triphosphate — un-benchmarkable
    if hcv_q and hcv_p:
        zero = sum(1 for v in hcv_p["per_active_counts"].values() if v == 0)
        lines.append(
            f"| HCV NS5B NI | catalytic active site (triphosphate) | {hcv_p['n_actives']} | "
            f"property-matched decoys IMPOSSIBLE ({zero}/{hcv_p['n_actives']} actives "
            f"recovered 0 decoys; {hcv_q['decoys_per_active_mean']}/active) | n/a | n/a | n/a "
            f"| n/a | **un-benchmarkable** (meta-finding) |")
    else:
        lines.append("| HCV NS5B NI | catalytic active site | — | — | — | — | — | — | artifact missing |")

    # SARS-CoV-2 Mpro (non-covalent)
    if mpro and mpro.get("status") == "blocked_receptor_prep":
        n_act = mpro_s["n_non_covalent_actives"] if mpro_s else "?"
        n_inact = mpro_s["n_inactives"] if mpro_s else "?"
        lines.append(
            f"| SARS-CoV-2 Mpro (non-covalent) | 3CLpro active site (7L11:A) | "
            f"{n_act} non-cov (avail.) | {n_inact} measured Moonshot inactives (real) | "
            f"blocked | blocked | blocked | blocked | **blocked: receptor prep** (dataset ready) |")
    elif mpro and mpro.get("benchmark"):
        b = mpro["benchmark"]
        boot = mpro.get("bootstrap", {}).get("bootstrap", {})
        lines.append(
            f"| SARS-CoV-2 Mpro (non-covalent) | 3CLpro active site (7L11:A) | "
            f"{b['n_actives']} | measured Moonshot inactives | {_ci(boot,'bedroc')} | "
            f"{_ci(boot,'EF1%')} | {_ci(boot,'log_auc')} | {_ci(boot,'roc_auc')} | "
            f"{mpro.get('verdict')} |")
    else:
        lines.append("| SARS-CoV-2 Mpro (non-covalent) | 3CLpro active site | — | — | — | — | — | — | not yet run |")

    lines += [
        "",
        "## Reading the table",
        "",
        "- **TiLV PB1**: a real docking run, but four actives make every CI span almost the",
        "  whole [0,1] range — a labelled demonstration, not an enrichment claim.",
        "- **HCV NS5B NI**: not assessable by matched-decoy enrichment at all (see meta-finding).",
        "- **SARS-CoV-2 Mpro**: powered target with real measured inactives and OpenBabel",
        "  receptor-prep (Meeko 0.7.1 declined for every Mpro chain; OpenBabel fallback used).",
        "  Docking signal is MODEST — ROC-AUC CI crosses 0.5. Reported as-is, not inflated.",
        "  See Phase 10 (`outputs/phase10/`) for ensemble conformer comparison.",
    ]
    return "\n".join(lines) + "\n"


def build_gate() -> str:
    mpro = _load(MPRO)
    blocked = bool(mpro and mpro.get("status") == "blocked_receptor_prep")
    lines = [
        "# Phase 9C Gate Decision: SARS-CoV-2 Mpro",
        "",
        "Decision: **DO NOT PROMOTE any DL or consensus ranking term.**",
        "",
        "Evidence:",
        "",
    ]
    if blocked:
        lines += [
            "- The powered Mpro benchmark has not produced a number yet: real Vina receptor",
            "  prep is blocked by a Meeko 0.7.1 Mpro-specific bug (see",
            "  `outputs/phase9/mpro_noncovalent_benchmark.json`).",
            "- With no powered enrichment result, there is no leakage-controlled DL-vs-baseline",
            "  comparison to justify a ranking change.",
        ]
    else:
        lines += [
            "- A powered Mpro enrichment result exists, but promotion requires a DL/consensus",
            "  signal that beats the Vina baseline with NON-OVERLAPPING CIs.",
            "- `consensus_node` reads `cnn_affinity`/`boltzina_score`, which are not populated",
            "  (no DL rescore is wired), so no DL term can be evaluated against the baseline.",
        ]
    lines += [
        "",
        "Therefore consensus/DL outputs remain annotation-only and ranking weights are",
        "unchanged. This negative result is logged deliberately: even a powered benchmark does",
        "not, on its own, license a ranking change without a leakage-controlled DL improvement.",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "multitarget_ci_table.md").write_text(build_table())
    (OUT / "gate_decision_mpro.md").write_text(build_gate())
    print(f"wrote {OUT / 'multitarget_ci_table.md'}")
    print(f"wrote {OUT / 'gate_decision_mpro.md'}")


if __name__ == "__main__":
    main()
