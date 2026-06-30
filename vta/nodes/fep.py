"""fep_node — opt-in FEP / ABFE validation seam (Phase 5).

Free-energy perturbation / absolute binding free energy is the most expensive
validation layer in the funnel. This node integrates it as a production hook without
pretending the fast path can run alchemical MD by default:

  * choose the top FEP_TOP_N MD-validated leads when available,
  * call an external ABFE runner if one is discoverable,
  * parse JSON results with delta_g / error fields,
  * otherwise label the candidates as skipped.

Runner contract:

    openfe-abfe <protein_pdb> <ligand_smiles> <output_dir>

The runner should print JSON to stdout, for example:

    {"delta_g": -8.4, "error": 0.7, "method": "OpenFE"}
"""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

from vta.state import VTAState

FEP_TOP_N = 3


def _fep_runner() -> str | None:
    from vta.toolconfig import find_tool
    return find_tool("openfe-abfe")


def _run_abfe(cmd: str, protein_pdb: str, smiles: str, out_dir: str) -> dict | None:
    """Run one ABFE job through the configured command and parse JSON stdout."""
    proc = subprocess.run(
        [cmd, protein_pdb, smiles, out_dir],
        capture_output=True, text=True, timeout=604800,
    )
    if proc.returncode != 0:
        return {"status": "failed", "reason": proc.stderr.strip() or "ABFE runner failed"}
    try:
        data = json.loads(proc.stdout)
    except ValueError:
        return {"status": "failed", "reason": "ABFE runner did not emit JSON"}
    delta_g = data.get("delta_g")
    if delta_g is None:
        return {"status": "failed", "reason": "ABFE JSON missing delta_g"}
    return {
        "status": "completed",
        "delta_g": round(float(delta_g), 3),
        "error": round(float(data["error"]), 3) if data.get("error") is not None else None,
        "method": data.get("method", "external ABFE runner"),
    }


def _candidate_pool(state: VTAState) -> list[dict]:
    return (
        state.get("md_validated_leads")
        or state.get("md_candidates")
        or state.get("lead_candidates")
        or []
    )


def fep_node(state: VTAState) -> VTAState:
    """Annotate top leads with FEP/ABFE results when a runner is available."""
    candidates = list(_candidate_pool(state))[:FEP_TOP_N]
    if not candidates:
        state["fep_results"] = {}
        state["fep_validated_leads"] = []
        state["audit_trail"].append("FEP: no candidates — phase skipped")
        return state

    runner = _fep_runner()
    if not runner:
        state["fep_results"] = {
            c["ligand"]: {"status": "skipped", "reason": "openfe-abfe not found"}
            for c in candidates
        }
        for c in candidates:
            c["fep_badge"] = "FEP not run"
        state["fep_validated_leads"] = candidates
        state["versions"]["fep"] = "skipped (no openfe-abfe)"
        state["audit_trail"].append(
            f"[skip] FEP: openfe-abfe not found — {len(candidates)} candidates labelled skipped")
        return state

    out_root = os.path.join("outputs", state["run_id"], "fep")
    os.makedirs(out_root, exist_ok=True)
    structures = state.get("structures") or {}
    results: dict[str, dict] = {}

    for c in candidates:
        ligand = c["ligand"]
        protein = c.get("protein")
        smiles = c.get("smiles")
        pdb_path = (structures.get(protein) or {}).get("pdb_path")
        if not (protein and smiles and pdb_path):
            results[ligand] = {
                "status": "skipped",
                "reason": "candidate missing protein PDB or SMILES",
            }
            c["fep_badge"] = "FEP skipped"
            continue

        run_dir = os.path.join(out_root, f"{ligand}_{protein}".replace(" ", "_"))
        os.makedirs(run_dir, exist_ok=True)
        result = _run_abfe(runner, pdb_path, smiles, run_dir)
        results[ligand] = result or {"status": "failed", "reason": "no result"}
        if results[ligand].get("status") == "completed":
            c["fep_delta_g"] = results[ligand]["delta_g"]
            c["fep_error"] = results[ligand].get("error")
            c["fep_badge"] = "FEP-computed"
        else:
            c["fep_badge"] = "FEP failed"

    ranked = sorted(candidates, key=lambda c: c.get("fep_delta_g", float("inf")))
    state["fep_results"] = results
    state["fep_validated_leads"] = ranked
    completed = sum(1 for r in results.values() if r.get("status") == "completed")
    state["versions"]["fep"] = "external ABFE runner"
    state["audit_trail"].append(
        f"FEP: {completed}/{len(candidates)} candidates completed via {Path(runner).name}")
    return state
