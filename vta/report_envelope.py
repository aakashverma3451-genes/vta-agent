"""The honesty envelope — the non-removable contract wrapped around every VTA-Agent output.

This is D0.5 of the deployment plan (the "report-level honesty contract"). Its job is to make
it **structurally impossible** to emit a ranking without the context a reader needs to not
misread it as an efficacy claim. `build_envelope()` is called unconditionally by
`report_node` (stored on the state) and by `render_report` (rendered as the first card), so
there is no code path — deferred run, empty run, or full run — that shows a lead without:

  1. the hypotheses-not-efficacy disclaimer (D0.5),
  2. the frozen benchmark the ranking's credibility rests on (grade + metric + 95% CI), pinned
     by content hash (D0.5 / D3 benchmark pinning),
  3. whether docking beat a **trivial 2-D-similarity baseline** by a *paired* test (D0.1/D0.3),
  4. the redocking **pose reliability** of the wired target (D0.2),
  5. the in-force ranking weights and *why* (ΔG-primary, LE demoted by the WI-6 gate — D0.4),
  6. the nucleotide/metal and modest-signal **scoring caveats** (Phase 9D),
  7. an **out-of-validated-domain** flag for any run target no frozen benchmark covers.

Every number is read from a committed artifact, never asserted inline: the frozen benchmark
(`outputs/phase9/locked_benchmark.json`), the paired-baseline + LE gate
(`outputs/phase11/baselines_and_gate.json`), and the redocking check
(`outputs/phase11/redock_validation.json`). Each read degrades to a *labelled* "not available"
rather than crashing — the honesty contract must never be the thing that takes a run down, but
it also must never silently omit itself. The live ranking weights are imported from
`vta.nodes.rank` so the envelope can never drift from what actually ranked.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from vta.nodes.rank import WEIGHTS as LIVE_RANK_WEIGHTS
from vta.state import VTAState

# Committed artifact locations (the ground truth the envelope cites). Overridable so tests are
# hermetic — point them at fixtures, or delete them to exercise the labelled-fallback path.
DEFAULT_BENCHMARK_PATH = Path("outputs/phase9/locked_benchmark.json")
DEFAULT_BASELINES_PATH = Path("outputs/phase11/baselines_and_gate.json")
DEFAULT_REDOCK_PATH = Path("outputs/phase11/redock_validation.json")

DISCLAIMER = (
    "Outputs are ranked, uncertainty-bearing HYPOTHESES for research triage — "
    "NOT clinical, efficacy, or safety claims. No result here has been experimentally "
    "confirmed. Every ranking rests on the frozen benchmark cited below; read its grade "
    "and 95% CI before acting on any lead."
)

# Map a run's protein name → the frozen-benchmark target, the baselines_and_gate benchmark key,
# and the redock_validation target id. Keyword-matched (case-insensitive substring) so that
# "MPRO", "Mpro_7L11", "SARS2-Mpro" all resolve. A run protein that matches none is flagged
# out-of-validated-domain rather than silently borrowing another target's grade.
_TARGET_MAP = (
    # (keywords, frozen_target_name, baselines_key, redock_id)
    (("mpro", "3clpro", "main protease"), "SARS-CoV-2 Mpro (non-covalent)",
     "Mpro_noncovalent", "Mpro_7L11"),
    (("ns5b", "hcv"), "HCV NS5B NI (triphosphate)", None, "HCV_NS5B_2XI3"),
    (("pb1", "tilv"), "TiLV PB1", "TiLV_PB1", "TiLV_8PSO"),
)


def _load_json(path: Path) -> Optional[Any]:
    """Read a committed artifact; return None (labelled upstream) on any failure."""
    try:
        return json.loads(Path(path).read_text())
    except Exception:
        return None


def _match_target(protein: str) -> Optional[tuple]:
    p = (protein or "").lower()
    for keywords, frozen, base_key, redock_id in _TARGET_MAP:
        if any(k in p for k in keywords):
            return frozen, base_key, redock_id
    return None


def _run_proteins(state: VTAState) -> List[str]:
    """Proteins that were actually ranked this run (leads first, then structures)."""
    seen: List[str] = []
    for lead in state.get("lead_candidates") or []:
        name = lead.get("protein")
        if name and name not in seen:
            seen.append(name)
    for name in (state.get("structures") or {}):
        if name not in seen:
            seen.append(name)
    return seen


def _frozen_target(benchmark: Optional[dict], name: str) -> Optional[dict]:
    for t in ((benchmark or {}).get("targets") or []):
        if t.get("target") == name:
            return t
    return None


def _baseline_verdict(baselines: Optional[Any], key: Optional[str]) -> Dict[str, Any]:
    """Did docking beat the trivial 2-D-similarity baseline by a PAIRED test? (D0.1/D0.3)."""
    if key is None:
        return {"available": False,
                "statement": "No powered matched-decoy benchmark exists for this target "
                             "(see scoring caveats); docking vs trivial-baseline is untested."}
    for entry in (baselines or []):
        if entry.get("benchmark") == key:
            beats = bool(entry.get("structure_based_enrichment_demonstrated"))
            return {
                "available": True,
                "beats_trivial_2d_baseline": beats,
                "statement": entry.get("wi3_plain_language")
                or ("docking beats 2D-similarity (paired test)" if beats
                    else "docking does NOT beat 2D-similarity (paired test)"),
            }
    return {"available": False,
            "statement": "Trivial-baseline comparison artifact not found — treat ranking as "
                         "unvalidated against a 2-D-similarity floor."}


def _pose_reliability(redock: Optional[dict], redock_id: Optional[str]) -> Dict[str, Any]:
    """Redocking / pose-reproduction status for the wired target (D0.2)."""
    for r in ((redock or {}).get("results") or []):
        if r.get("target") == redock_id:
            reliable = r.get("pose_reliable")
            rmsd = r.get("rmsd")
            if reliable is True:
                status = "pose-reliable"
            elif reliable is False:
                status = "pose-UNRELIABLE"
            else:
                status = "pose reliability UNKNOWN (RMSD uncomputed)"
            return {"available": True, "status": status, "rmsd_A": rmsd,
                    "threshold_A": (redock or {}).get("rmsd_flag_threshold_A"),
                    "detail": r.get("status")}
    return {"available": False,
            "status": "pose reliability UNKNOWN (no redocking record for this target)"}


def _scoring_caveats(frozen: Optional[dict], names: List[str]) -> List[str]:
    caveats: List[str] = []
    # Always carry the class-level nucleotide/metal caveat established in Phase 9D — it warns
    # a reader that nucleotide-analog antivirals are scored on their parent, not active form.
    if any(_match_target(n) and _match_target(n)[0] == "HCV NS5B NI (triphosphate)"
           for n in names):
        caveats.append(
            "NUCLEOTIDE/METAL caveat: nucleotide-analog antivirals act as charged "
            "triphosphates with no property-matched decoys; docking here scores the parent "
            "species, not the active form (Phase 9D). A ranking for this class is NOT a "
            "validated enrichment.")
    # Pull the frozen benchmark's own honest caveats (modest signal, underpowered targets).
    for c in ((frozen or {}).get("caveats") or []):
        caveats.append(c)
    return caveats


def build_envelope(
    state: VTAState,
    *,
    benchmark_path: Path = DEFAULT_BENCHMARK_PATH,
    baselines_path: Path = DEFAULT_BASELINES_PATH,
    redock_path: Path = DEFAULT_REDOCK_PATH,
) -> Dict[str, Any]:
    """Build the honesty envelope for a run. Pure w.r.t. the state; reads committed artifacts.

    Always returns a dict carrying the disclaimer and the live ranking weights, even when every
    artifact is missing (each missing piece is labelled, never fabricated or omitted silently).
    """
    benchmark = _load_json(benchmark_path)
    baselines = _load_json(baselines_path)
    redock = _load_json(redock_path)

    if benchmark:
        pin = {
            "available": True,
            "artifact": benchmark.get("artifact"),
            "frozen_at": benchmark.get("frozen_at"),
            "content_hash": benchmark.get("content_hash"),
            "gate_decision": benchmark.get("gate_decision"),
        }
    else:
        pin = {"available": False,
               "note": "frozen benchmark artifact not found — ranking has no pinned "
                       "validation basis; treat as exploratory."}

    ranking_basis: List[Dict[str, Any]] = []
    out_of_domain: List[str] = []
    proteins = _run_proteins(state)
    for name in proteins:
        matched = _match_target(name)
        if not matched:
            out_of_domain.append(name)
            continue
        frozen_name, base_key, redock_id = matched
        ft = _frozen_target(benchmark, frozen_name)
        ranking_basis.append({
            "run_protein": name,
            "benchmark_target": frozen_name,
            "in_validated_domain": True,
            "grade": (ft or {}).get("grade"),
            "metrics_ci": (ft or {}).get("metrics_ci"),
            "trivial_baseline": _baseline_verdict(baselines, base_key),
            "pose_reliability": _pose_reliability(redock, redock_id),
        })

    envelope: Dict[str, Any] = {
        "disclaimer": DISCLAIMER,
        "benchmark_pin": pin,
        "ranking_weights": dict(LIVE_RANK_WEIGHTS),
        "ranking_note": (
            "ΔG-primary ranking. Ligand efficiency and conservation are REPORTED annotations "
            "(weight 0), not ranking terms — LE was demoted by the Phase 11 WI-6 paired-test "
            "gate (the LE-led composite did not beat ΔG-only on the powered Mpro benchmark)."),
        "ranking_basis": ranking_basis,
        "out_of_validated_domain": out_of_domain,
        "scoring_caveats": _scoring_caveats(benchmark, proteins),
        "provenance": {
            "frozen_benchmark": str(benchmark_path) if benchmark else "MISSING",
            "baselines_and_gate": str(baselines_path) if baselines else "MISSING",
            "redock_validation": str(redock_path) if redock else "MISSING",
        },
    }
    return envelope
