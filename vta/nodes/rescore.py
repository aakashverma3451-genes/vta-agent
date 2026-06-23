"""rescore_node — deep-learning re-scoring of docked poses (Phase 3b seam).

Re-scores each Vina pose with GNINA's CNN scoring function, surfacing two new
signals per docking record:

    cnn_score      pose-quality probability (0-1, higher = more pose-like)
    cnn_affinity   predicted binding affinity (pKd units, higher = tighter)

Why this matters. Scoring is the pipeline's most sensitive layer — we already proved
it once (ΔG-led → LE-led flipped the validation gate). Vina's empirical ΔG correlates
only modestly with experiment and carries a molecular-size bias we had to correct in
ranking; a learned scorer like GNINA's CNN correlates better and is the highest-ceiling
accuracy gain in the integration plan.

ANNOTATION-ONLY, DELIBERATELY. This node writes new `cnn_*` fields; it does NOT feed
the composite in rank.py yet. Blending a DL term into the score re-opens the
validation-gate calibration (the same discipline that caught the ΔG size bias), and
GNINA's CNN really wants a CUDA GPU this arm64 box doesn't have. So: build the seam
now (auto-detect, score real poses where possible, store the fields), blend the term
and re-run the gate on a GPU host later.

Auto-detects the `gnina` binary and degrades to a labelled skip (Vina ΔG retained)
when absent — same don't-crash discipline as structure/pockets/docking/admet.
"""
from __future__ import annotations

import os
import re
import subprocess

from vta.state import VTAState

# gnina --score_only prints these; higher is better for both.
_CNN_SCORE = re.compile(r"CNNscore:\s*(-?\d+\.?\d*)")
_CNN_AFFINITY = re.compile(r"CNNaffinity:\s*(-?\d+\.?\d*)")


def _gnina_bin() -> str | None:
    from vta.toolconfig import find_tool
    return find_tool("gnina", "GNINA_BIN")


def _run_gnina(gnina: str, receptor: str, pose: str) -> tuple[float, float] | None:
    """Re-score one docked pose; return (cnn_score, cnn_affinity) or None on failure.

    `--score_only` evaluates the pose as-is (no re-docking), so this is a cheap
    rescore of the Vina output, not a second search.
    """
    p = subprocess.run(
        [gnina, "--receptor", receptor, "--ligand", pose, "--score_only"],
        capture_output=True, text=True, timeout=600,
    )
    s = _CNN_SCORE.search(p.stdout)
    a = _CNN_AFFINITY.search(p.stdout)
    if not s or not a:
        return None
    return float(s.group(1)), float(a.group(1))


def rescore_node(state: VTAState) -> VTAState:
    """Re-score Vina poses with GNINA's CNN when available, else labelled skip."""
    rows = state.get("docking_results") or []
    if not rows:
        return state

    gnina = _gnina_bin()
    if not gnina:
        state["audit_trail"].append(
            "[skip] DL-rescore: gnina not installed (Vina ΔG retained)")
        state["versions"]["rescore"] = "skipped"
        return state

    rescored = 0
    for r in rows:
        # Only the real Vina path saves pose/receptor PDBQTs; mock records have none.
        pose, receptor = r.get("pose_path"), r.get("receptor_path")
        if not (pose and receptor and os.path.exists(pose) and os.path.exists(receptor)):
            continue
        try:
            out = _run_gnina(gnina, receptor, pose)
        except Exception as e:
            state["audit_trail"].append(
                f"DL-rescore: {r.get('ligand')} error ({type(e).__name__}: {e}); skipped")
            continue
        if out is None:
            continue
        r["cnn_score"], r["cnn_affinity"] = round(out[0], 3), round(out[1], 3)
        rescored += 1

    if rescored:
        state["versions"]["rescore"] = "GNINA (CNN)"
        best = max((r for r in rows if "cnn_affinity" in r),
                   key=lambda r: r["cnn_affinity"])
        state["audit_trail"].append(
            f"DL-rescore: GNINA CNN re-scored {rescored} poses (annotation-only, "
            f"ranking unchanged); best CNNaffinity {best['cnn_affinity']} "
            f"({best.get('ligand')})")
    else:
        # gnina present but nothing scoreable (e.g. mock docking, no saved poses).
        state["audit_trail"].append(
            "[skip] DL-rescore: gnina present but no docked poses to score")
        state["versions"]["rescore"] = "skipped (no poses)"
    return state
