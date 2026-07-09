"""run_triage — the honest triage a chat turn runs against VTA-Agent.

A chat message names a target (a genome/FASTA needs the full folding+docking pipeline and is
out of scope for an interactive turn — we say so). We run the REAL, offline, fast reasoning
layer — the same committed nodes the graph uses — and ground the verdict in the frozen
benchmark, so nothing here is fabricated:

  dossier_node        → structure provenance, pocket, benchmarkability (from locked_benchmark)
  triage_router_node  → full_dock | annotate_only | defer | refuse (+ rationale)
  build_envelope      → the non-removable honesty envelope: pinned benchmark grade + 95% CI,
                        docking-vs-2D-baseline paired verdict, pose reliability, caveats

The reply always carries the disclaimer (ranked hypotheses, not efficacy claims) and, for the
powered Mpro target, the powered fair-arena cliff result if the artifact is present. No docking
runs and no network is touched — an interactive turn stays honest and fast.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from vta.nodes.dossier import dossier_node
from vta.nodes.triage import triage_router_node
from vta.report_envelope import build_envelope
from vta.state import new_state

# Curated experimental structures for the benchmarked targets (mirrors the pipeline's wiring;
# domain facts, not fabrication). An unnamed target folds to a predicted stub → out-of-domain.
_KNOWN = {
    ("mpro", "3clpro", "main protease", "sars"): ("MPRO", {"method": "experimental",
                                                            "source": "7L11:A", "mean_plddt": None}),
    ("ns5b", "hcv"): ("NS5B", {"method": "experimental", "source": "2XI3:A", "mean_plddt": None}),
    ("pb1", "tilv", "polymerase", "rdrp"): ("PB1", {"method": "experimental",
                                                    "source": "8PSO:B", "mean_plddt": None}),
}
_CLIFF = Path("outputs/phaseS/activity_cliffs.json")
_RFCLIFF = Path("outputs/phaseS/rfscore_cliff.json")


def _resolve_target(query: str) -> Tuple[str, Dict[str, Any], bool]:
    """Map a free-text query to (protein_name, structure_record, is_known)."""
    q = (query or "").lower()
    for keys, (name, rec) in _KNOWN.items():
        if any(k in q for k in keys):
            return name, dict(rec), True
    # unknown → treat a short token as a protein name with a predicted-structure stub
    token = "".join(c for c in (query or "").strip().split("\n")[0][:32] if c.isalnum() or c in "-_") or "TARGET"
    return token.upper(), {"method": "esmfold", "mean_plddt": None}, False


def _looks_like_sequence(query: str) -> bool:
    s = "".join((query or "").split())
    if len(s) < 40:
        return False
    aa = sum(c.upper() in "ACDEFGHIKLMNPQRSTVWY" for c in s)
    return aa / max(1, len(s)) > 0.9 or (query or "").lstrip().startswith(">")


def _cliff_line() -> Optional[str]:
    try:
        c = json.loads(_CLIFF.read_text())["results"]
        vina, twod = c["vina_accuracy"], c["twod_knn_accuracy"]
        rf = None
        if _RFCLIFF.exists():
            rf = json.loads(_RFCLIFF.read_text())["accuracies"]["rfscore"]
        line = (f"**Fair-arena test (activity cliffs, {c['n_pairs']} pairs, chance 0.50):** "
                f"Vina {vina['median']} {vina['ci95']} vs 2-D-kNN {twod['median']} {twod['ci95']}")
        if rf:
            line += f"; learned RF-Score {rf['median']} {rf['ci95']}"
        return line + " — docking is below chance; the trivial ligand baseline wins."
    except Exception:
        return None


def run_triage(query: str) -> Dict[str, Any]:
    """Return {target, decision, steps[], markdown} for a chat query. Offline, no docking."""
    if _looks_like_sequence(query):
        return {
            "target": None, "decision": "defer",
            "steps": ["Input looks like a raw sequence/FASTA."],
            "markdown": (
                "**A raw genome/sequence needs the full folding + docking pipeline** (structure "
                "prediction, pocket detection, real AutoDock Vina) — that is not run in an "
                "interactive chat turn. Name a benchmarked target instead (e.g. *SARS-CoV-2 "
                "Mpro*, *HCV NS5B*, *TiLV PB1*) for the honest triage verdict, or run the full "
                "pipeline offline via `build_app`.\n\n"
                "_Outputs are ranked, uncertainty-bearing hypotheses — never clinical, efficacy, "
                "or safety claims._"),
        }

    name, rec, known = _resolve_target(query)
    st = new_state(f"chat-{name}", f"chat-{name}")
    st["structures"] = {name: rec}
    st = dossier_node(st)
    st = triage_router_node(st)
    decision = st["triage_decision"][name]["decision"]
    rationale = "; ".join(st["triage_decision"][name].get("rationale") or [])
    env = build_envelope(st)

    lines: List[str] = [f"## VTA-Agent — triage for **{name}**", "",
                        f"**Routing decision:** `{decision}`" + (f" — {rationale}" if rationale else "")]

    basis = next((b for b in env.get("ranking_basis") or [] if b.get("run_protein") == name), None)
    if basis:
        tb = basis.get("trivial_baseline") or {}
        pr = basis.get("pose_reliability") or {}
        mci = basis.get("metrics_ci") or {}
        lines += ["", f"**Benchmark grade:** {basis.get('grade')}"]
        if tb.get("statement"):
            beats = tb.get("beats_trivial_2d_baseline")
            mark = "✗ does NOT beat" if beats is False else ("✓ beats" if beats else "untested vs")
            lines.append(f"**Docking vs trivial 2-D baseline:** {mark} — {tb['statement']}")
        if mci:
            b, r = mci.get("BEDROC") or {}, mci.get("ROC_AUC") or {}
            if b:
                lines.append(f"**BEDROC(α20):** {b.get('median')} {b.get('ci95')}")
            if r:
                lines.append(f"**ROC-AUC:** {r.get('median')} {r.get('ci95')}")
        lines.append(f"**Pose reliability:** {pr.get('status')}"
                     + (f" ({pr.get('rmsd_A')} Å)" if pr.get('rmsd_A') is not None else ""))
        if name == "MPRO":
            cl = _cliff_line()
            if cl:
                lines += ["", cl]
    elif env.get("out_of_validated_domain"):
        lines += ["", "⚠ **Out-of-validated-domain** — no frozen benchmark covers this target; "
                  "any ranking would be exploratory and is not backed by a validation grade."]

    caveats = env.get("scoring_caveats") or []
    if caveats:
        lines += ["", "**Scoring caveats:**"] + [f"- {c}" for c in caveats[:3]]

    pin = env.get("benchmark_pin") or {}
    lines += ["", "---", f"_{env.get('disclaimer')}_"]
    if pin.get("available"):
        lines.append(f"_Pinned benchmark: {pin.get('artifact')} · frozen {pin.get('frozen_at')} · "
                     f"hash {pin.get('content_hash')} · gate: {pin.get('gate_decision')}_")

    return {"target": name, "decision": decision,
            "steps": list(st.get("audit_trail") or []),
            "markdown": "\n".join(lines)}
