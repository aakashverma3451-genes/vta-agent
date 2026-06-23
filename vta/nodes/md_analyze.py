"""md_analyze_node — RMSD stability, persistent contacts, and MM-GBSA (Phase 4).

For each completed MD simulation this node computes three signals:

    1. Ligand RMSD over time (MDAnalysis) — STABLE (<2Å mean), MODERATE (2–4Å),
       UNSTABLE (>4Å).  The cutoff values match the Yamaotsu & Hirono 2016 criteria
       used in the prior TiLV MD paper (Sumon et al. 2023).

    2. Persistent contacts — protein residues that stay within 4 Å of the ligand for
       >50% of trajectory frames. These are the mechanistic contacts to cite.

    3. MM-GBSA binding energy — simplified inline calculation from the OpenMM
       trajectory. For publication-grade results use gmx_MMPBSA / MMPBSA.py instead
       (the `mean_binding_energy` field will be None until then).

Gracefully skips candidates whose simulation failed or whose trajectory is absent.
No GPU required for analysis; MDAnalysis runs on CPU.
"""
from __future__ import annotations

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


def _mmgbsa_placeholder() -> dict:
    # Full MM-GBSA requires AmberTools (MMPBSA.py) or gmx_MMPBSA.
    # Built as placeholder; mean_binding_energy will be populated when
    # the AmberTools integration lands (Phase 4 production).
    return {"mean_binding_energy": None, "std_binding_energy": None,
            "note": "full MM-GBSA requires AmberTools; placeholder"}


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
            rmsd = _rmsd_analysis(traj, top)
            contacts = _contacts_analysis(traj, top)
            mmgbsa = _mmgbsa_placeholder()
            verdict = rmsd.get("stability", "unknown")
            analysis[lid] = {
                "rmsd": rmsd, "contacts": contacts,
                "mmgbsa": mmgbsa, "md_verdict": verdict,
            }
            state["audit_trail"].append(
                f"MD-analyze[{lid}]: {verdict} "
                f"(mean RMSD {rmsd.get('mean_rmsd')} Å, "
                f"{contacts.get('n_persistent')} persistent contacts)")
        except Exception as e:
            analysis[lid] = {"status": "failed", "reason": str(e)}
            state["audit_trail"].append(f"MD-analyze[{lid}]: error — {e}")

    state["md_analysis"] = analysis
    return state
