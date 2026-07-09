"""verification_node (Phase R / R4) — the hard gate every claim passes before the report.

HemaGuide's retrieval-first, tool-constrained design cut hallucination to 0.3%; the decisive
verification for a physics pipeline is deterministic physics + statistics, not LLM self-critique.
This node sits on the single edge into reporting and writes `verification_verdict`; the report
obeys it (a docking ranking that fails the gate is downgraded to a labelled annotation).

Three deterministic gates decide whether a *docking-enrichment* claim may stand:

  1. redock-RMSD gate  — the wired target's co-crystal ligand must redock < 2.0 Å (Trott & Olson
                         2010; Quiroga & Villarreal 2016). Read from redock_validation.json.
  2. paired-baseline   — Vina must SIGNIFICANTLY beat the trivial 2-D-similarity baseline by the
                         paired-bootstrap test (Holm). Read `structure_based_enrichment_demonstrated`
                         from baselines_and_gate.json. This is the D0.3 gate promoted to a HARD
                         pre-report gate.
  3. applicability     — the target must be in-domain for the scorer (triage full_dock). Targets
                         routed annotate_only/defer/refuse are not making a docking claim.

Verdict (run-level, worst-case across targets that carry a docking-based lead):
  pass       — every docking claim clears gates 1–3.
  downgrade  — a docking claim fails the paired-baseline (or pose) gate → present as annotation.
  defer/refuse — the router deferred/refused (no docking claim to verify).
  not_applicable — a deferred run with no leads.

An LLM chain-of-verification layer over the narrative may be added later; per design it may
NEVER override these gates (Zheng et al. 2023 documents LLM-judge biases).
"""
from __future__ import annotations

from typing import Any, Dict, List

from vta.report_envelope import (
    DEFAULT_BASELINES_PATH,
    DEFAULT_REDOCK_PATH,
    _baseline_verdict,
    _load_json,
    _match_target,
    _pose_reliability,
)
from vta.state import VTAState


def _lead_proteins(state: VTAState) -> List[str]:
    seen: List[str] = []
    for lead in state.get("lead_candidates") or []:
        p = lead.get("protein")
        if p and p not in seen:
            seen.append(p)
    return seen


def build_verdict(state: VTAState, *, baselines_path=DEFAULT_BASELINES_PATH,
                  redock_path=DEFAULT_REDOCK_PATH) -> Dict[str, Any]:
    """Pure: compute the verification verdict from the state + committed gate artifacts."""
    baselines = _load_json(baselines_path)
    redock = _load_json(redock_path)
    triage = state.get("triage_decision") or {}
    lead_proteins = _lead_proteins(state)

    per_target: Dict[str, Any] = {}
    worst = "pass"
    order = {"pass": 0, "not_applicable": 0, "defer": 1, "refuse": 2, "downgrade": 3}

    # Targets making a docking-enrichment claim: those with docking-based leads AND full_dock.
    for protein in lead_proteins:
        decision = (triage.get(protein) or {}).get("decision", "full_dock")
        matched = _match_target(protein)
        base = _baseline_verdict(baselines, matched[1] if matched else None)
        pose = _pose_reliability(redock, matched[2] if matched else None)
        beats = base.get("beats_trivial_2d_baseline")
        pose_ok = pose.get("status") == "pose-reliable"

        if decision != "full_dock":
            tv = "downgrade"          # routed away from docking but leads exist → annotate
            why = f"triage routed {protein} to {decision}"
        elif beats is True and pose_ok:
            tv = "pass"
            why = "beats 2D baseline (paired) and pose-reliable"
        elif beats is False:
            tv = "downgrade"
            why = "docking does NOT beat the 2D-similarity baseline (paired test)"
        elif not pose_ok:
            tv = "downgrade"
            why = f"pose reliability not established ({pose.get('status')})"
        else:
            tv = "downgrade"
            why = "docking-vs-2D baseline untested for this target — enrichment unverifiable"
        per_target[protein] = {"verdict": tv, "decision": decision, "beats_2d_baseline": beats,
                               "pose": pose.get("status"), "rationale": why}
        if order[tv] > order[worst]:
            worst = tv

    # Targets routed defer/refuse (no docking claim) — record for completeness.
    for protein, d in triage.items():
        dec = (d or {}).get("decision")
        if dec in {"defer", "refuse"} and protein not in per_target:
            per_target[protein] = {"verdict": dec, "decision": dec,
                                   "rationale": f"router {dec} — no docking claim"}
            if order[dec] > order[worst]:
                worst = dec

    if not lead_proteins and not per_target:
        worst = "not_applicable"

    gates = {
        "redock_rmsd": {"source": str(redock_path), "available": bool(redock)},
        "paired_baseline": {"source": str(baselines_path), "available": bool(baselines),
                            "rule": "Vina must beat 2D-similarity by paired bootstrap (Holm)"},
        "applicability_domain": {"rule": "target must be full_dock (in-domain) to make a "
                                         "docking claim"},
    }
    return {"verdict": worst, "gates": gates, "per_target": per_target,
            "note": ("A docking-enrichment claim is emitted only on `pass`; otherwise the "
                     "ranking is downgraded to a labelled ligand-based annotation.")}


def verification_node(state: VTAState) -> VTAState:
    verdict = build_verdict(state)
    state["verification_verdict"] = verdict
    downgraded = [p for p, v in verdict["per_target"].items() if v["verdict"] == "downgrade"]
    state["audit_trail"].append(
        f"Verification: run verdict = {verdict['verdict'].upper()}"
        + (f"; downgraded docking claims for {', '.join(sorted(downgraded))}" if downgraded else ""))
    state["versions"]["verification"] = "physics+stats hard gate v1 (R4)"
    return state
