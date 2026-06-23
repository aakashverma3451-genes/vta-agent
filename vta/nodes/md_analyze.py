"""md_analyze_node — RMSD stability, persistent contacts, and MM-GBSA (Phase 4).

For each completed MD simulation this node computes three signals:

    1. Ligand RMSD over time (MDAnalysis) — STABLE (<2Å mean), MODERATE (2–4Å),
       UNSTABLE (>4Å).  The cutoff values match the Yamaotsu & Hirono 2016 criteria
       used in the prior TiLV MD paper (Sumon et al. 2023).

    2. Persistent contacts — protein residues that stay within 4 Å of the ligand for
       >50% of trajectory frames. These are the mechanistic contacts to cite.

    3. MM-GBSA binding free energy (ΔTOTAL) — computed by gmx_MMPBSA (preferred) or
       AmberTools' MMPBSA.py over the trajectory, parsed from the standard
       FINAL_RESULTS_MMPBSA.dat. This is the quantitative binding number reviewers
       expect alongside RMSD stability (cf. Sumon et al. 2023). Auto-detects the
       engine and degrades to a labelled skip (mean_binding_energy=None) when it —
       or the Amber topology it needs — is absent, so the MD report falls back to
       stability-only rather than crashing or fabricating a number.

ANNOTATION-ONLY: the energy is surfaced in analysis/report but NOT folded into
md_rerank's score, which still ranks on the RMSD verdict. Blending ΔG into ranking
requires a validation-gate recalibration (same discipline as rescore/boltzina), so
it's a deliberate later step, not a silent one here.

Gracefully skips candidates whose simulation failed or whose trajectory is absent.
No GPU required for analysis; MDAnalysis runs on CPU.
"""
from __future__ import annotations

import os
import re
import subprocess

from vta.state import VTAState

_STABLE_CUTOFF = 2.0    # Å — mean RMSD threshold for STABLE verdict
_MODERATE_CUTOFF = 4.0  # Å — mean RMSD threshold for MODERATE (else UNSTABLE)
_CONTACT_RADIUS = 4.0   # Å — protein–ligand contact definition
_CONTACT_OCC = 0.5      # fraction of frames required for "persistent"
_SAMPLE_EVERY = 10      # analyse every N-th trajectory frame (speed vs accuracy)


def _check_mdanalysis() -> bool:
    try:
        import MDAnalysis  # noqa: F401
        return True
    except ImportError:
        return False


def _rmsd_analysis(traj: str, topology: str) -> dict:
    import numpy as np
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms

    u = mda.Universe(topology, traj)
    ligand = u.select_atoms("resname LIG MOL UNK")
    if not ligand:
        return {"error": "ligand selection empty — check resname in PDB"}
    backbone = u.select_atoms("protein and name CA")
    R = rms.RMSD(ligand, ligand, ref_frame=0, superposition_group=backbone)
    R.run()
    vals = R.results.rmsd[:, 2]
    mean_r = float(np.mean(vals))
    return {
        "mean_rmsd": round(mean_r, 2),
        "max_rmsd": round(float(np.max(vals)), 2),
        "std_rmsd": round(float(np.std(vals)), 2),
        "final_rmsd": round(float(vals[-1]), 2),
        "stability": ("STABLE" if mean_r < _STABLE_CUTOFF
                      else "MODERATE" if mean_r < _MODERATE_CUTOFF
                      else "UNSTABLE"),
        "rmsd_timeseries": [round(float(v), 3) for v in vals],
    }


def _contacts_analysis(traj: str, topology: str) -> dict:
    import MDAnalysis as mda

    u = mda.Universe(topology, traj)
    ligand = u.select_atoms("resname LIG MOL UNK")
    total_sampled = 0
    counts: dict[str, int] = {}
    for _ in u.trajectory[::_SAMPLE_EVERY]:
        total_sampled += 1
        nearby = u.select_atoms(
            f"protein and around {_CONTACT_RADIUS} group lig", lig=ligand)
        for res in nearby.residues:
            key = f"{res.resname}{res.resid}"
            counts[key] = counts.get(key, 0) + 1

    if not total_sampled:
        return {"persistent_contacts": {}, "n_persistent": 0, "top_contacts": []}

    persistent = {k: round(v / total_sampled, 2)
                  for k, v in counts.items()
                  if v / total_sampled >= _CONTACT_OCC}
    top = sorted(persistent.items(), key=lambda x: -x[1])[:10]
    return {"persistent_contacts": persistent,
            "n_persistent": len(persistent),
            "top_contacts": top}


# Standard single-trajectory GB input (GB-OBC2, igb=2, physiological salt).
_MMGBSA_INPUT = (
    "MM-GBSA over the production trajectory\n"
    "&general\n  startframe=1, interval=10, verbose=2,\n/\n"
    "&gb\n  igb=2, saltcon=0.150,\n/\n"
)


def _mmgbsa_bin() -> str | None:
    """Find the MM-GBSA engine: gmx_MMPBSA (preferred) or AmberTools' MMPBSA.py."""
    from vta.toolconfig import find_tool
    return (find_tool("gmx_MMPBSA", "GMX_MMPBSA_BIN")
            or find_tool("MMPBSA.py", "MMPBSA_BIN"))


def _parse_mmgbsa_dat(text: str) -> dict | None:
    """Pull ΔTOTAL (binding free energy) out of a FINAL_RESULTS_MMPBSA.dat.

    Both MMPBSA.py and gmx_MMPBSA write, in the
    'Differences (Complex - Receptor - Ligand)' block, a line like:
        DELTA TOTAL          -28.45        3.12        0.10
    (Average, Std. Dev., Std. Err.). gmx_MMPBSA may render it as 'ΔTOTAL'.
    Returns {mean_binding_energy, std_binding_energy} in kcal/mol, or None.
    """
    for line in text.splitlines():
        norm = line.strip().replace("Δ", "DELTA ").upper()
        if norm.startswith("DELTA TOTAL"):
            nums = re.findall(r"-?\d+\.\d+", line)
            if not nums:
                return None
            return {
                "mean_binding_energy": round(float(nums[0]), 2),
                "std_binding_energy": round(float(nums[1]), 2) if len(nums) > 1 else None,
            }
    return None


def _run_mmgbsa_tool(binpath: str, parm: str, traj: str, run_dir: str) -> str | None:
    """Invoke the MM-GBSA engine; return FINAL_RESULTS_MMPBSA.dat text or None.

    Isolated as its own seam so tests can stub the external run without a binary.
    """
    os.makedirs(run_dir, exist_ok=True)
    in_path = os.path.join(run_dir, "mmgbsa.in")
    out_path = os.path.join(run_dir, "FINAL_RESULTS_MMPBSA.dat")
    with open(in_path, "w") as fh:
        fh.write(_MMGBSA_INPUT)
    proc = subprocess.run(
        [binpath, "-O", "-i", in_path, "-cp", parm, "-y", traj, "-o", out_path],
        capture_output=True, text=True, timeout=7200, cwd=run_dir,
    )
    if proc.returncode != 0 or not os.path.exists(out_path):
        return None
    with open(out_path) as fh:
        return fh.read()


def _compute_mmgbsa(result: dict, run_dir: str) -> dict:
    """Real MM-GBSA binding free energy with honest, labelled graceful fallback.

    Always returns a dict carrying `mean_binding_energy` (None when unavailable) so
    the report shape is stable. MM-GBSA needs an Amber topology (prmtop); the current
    OpenMM md_simulate doesn't emit one yet, so when the engine is present but no
    `parm`/`prmtop` is in the result we say so explicitly instead of faking a value.
    """
    binpath = _mmgbsa_bin()
    if not binpath:
        return {"mean_binding_energy": None, "std_binding_energy": None,
                "note": "[skip] MM-GBSA: gmx_MMPBSA / AmberTools (MMPBSA.py) not installed"}

    parm = result.get("parm") or result.get("prmtop")
    traj = result.get("trajectory")
    if not (parm and traj):
        return {"mean_binding_energy": None, "std_binding_energy": None,
                "note": "MM-GBSA engine present but no Amber topology (prmtop) in "
                        "md_results — md_simulate must emit one (ParmEd) for production"}

    try:
        dat = _run_mmgbsa_tool(binpath, parm, traj, run_dir)
    except Exception as e:
        return {"mean_binding_energy": None, "std_binding_energy": None,
                "note": f"MM-GBSA run failed: {e}"}
    if not dat:
        return {"mean_binding_energy": None, "std_binding_energy": None,
                "note": "MM-GBSA run produced no parseable results file"}

    parsed = _parse_mmgbsa_dat(dat)
    if not parsed:
        return {"mean_binding_energy": None, "std_binding_energy": None,
                "note": "MM-GBSA results file present but ΔTOTAL not found"}
    return {**parsed, "method": os.path.basename(binpath),
            "note": "MM-GBSA (GB-OBC2, igb=2) binding free energy"}


def md_analyze_node(state: VTAState) -> VTAState:
    """Analyse completed MD trajectories; skip failed/absent ones gracefully."""
    md_results = state.get("md_results") or {}
    if not md_results:
        state["md_analysis"] = {}
        return state

    if not _check_mdanalysis():
        state["md_analysis"] = {
            lid: {"status": "skipped", "reason": "MDAnalysis not installed"}
            for lid in md_results
        }
        state["audit_trail"].append(
            "[skip] MD-analyze: MDAnalysis not installed — install mdanalysis")
        return state

    analysis: dict = {}
    for lid, result in md_results.items():
        if result.get("status") != "completed":
            analysis[lid] = {"status": "skipped",
                             "reason": result.get("reason", "simulation not completed")}
            continue
        try:
            traj = result["trajectory"]
            top = result["final_pdb"]
            run_dir = os.path.dirname(traj) or "."
            rmsd = _rmsd_analysis(traj, top)
            contacts = _contacts_analysis(traj, top)
            mmgbsa = _compute_mmgbsa(result, run_dir)
            verdict = rmsd.get("stability", "unknown")
            analysis[lid] = {
                "rmsd": rmsd, "contacts": contacts,
                "mmgbsa": mmgbsa, "md_verdict": verdict,
            }
            dG = mmgbsa.get("mean_binding_energy")
            energy_str = f", MM-GBSA {dG} kcal/mol" if dG is not None else ""
            state["audit_trail"].append(
                f"MD-analyze[{lid}]: {verdict} "
                f"(mean RMSD {rmsd.get('mean_rmsd')} Å, "
                f"{contacts.get('n_persistent')} persistent contacts{energy_str})")
        except Exception as e:
            analysis[lid] = {"status": "failed", "reason": str(e)}
            state["audit_trail"].append(f"MD-analyze[{lid}]: error — {e}")

    state["md_analysis"] = analysis
    return state
