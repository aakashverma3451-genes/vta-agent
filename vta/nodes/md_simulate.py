"""md_simulate_node — run 100 ns OpenMM MD on each selected lead (Phase 4).

For each md_candidate this node:
  1. Converts the docked SMILES → 3D SDF (RDKit ETKDG)
  2. Builds a solvated protein-ligand system (OpenMMDL / OpenFF 2.0)
  3. Runs energy minimisation → 1 ns NVT equilibration → 100 ns NPT production
  4. Saves DCD trajectory + CSV energy log

Level-1 smoke mode (CPU / short run): set MD_STEPS=5_000_000 (10 ns) and
EQUILIBRATION_STEPS=100_000 (0.2 ns) via env vars MD_STEPS / MD_EQUIL_STEPS
to get a fast CPU demo. Production requires a CUDA GPU.

Auto-detects OpenMM and degrades to a labelled skip when absent — the same
don't-crash discipline as every other real node. CI never loads OpenMM.
"""
from __future__ import annotations

import os
from pathlib import Path

from vta.state import VTAState

# Env-var overrides let the smoke-test use 10 ns without code changes.
_MD_STEPS = int(os.environ.get("MD_STEPS", 50_000_000))       # 100 ns at 2 fs
_EQUIL_STEPS = int(os.environ.get("MD_EQUIL_STEPS", 500_000))  # 1 ns
_REPORT_INTERVAL = int(os.environ.get("MD_REPORT_INTERVAL", 50_000))  # every 100 ps
_TEMPERATURE = 300      # K
_SALT_CONC = 0.15       # M NaCl
_DT_FS = 0.002          # ps per step


def _check_openmm() -> bool:
    try:
        import openmm  # noqa: F401
        return True
    except ImportError:
        return False


def _check_parmed() -> bool:
    try:
        import parmed  # noqa: F401
        return True
    except ImportError:
        return False


def _save_amber_topology(system, topology, positions, out_dir: str, run_id: str) -> str | None:
    """Write an Amber complex topology (prmtop + inpcrd) via ParmEd, for MM-GBSA.

    md_analyze's MM-GBSA reads `result["parm"]`, but OpenMM doesn't speak Amber natively,
    so production runs had no topology to feed it. ParmEd translates the PARAMETERISED
    OpenMM System (forces set — `_run_simulation` even added the barostat to it) into
    Amber prmtop/inpcrd. Returns the prmtop path, or None on any failure — the caller
    then omits `parm` and md_analyze prints its own honest "no topology" note (so we
    don't duplicate it here).
    """
    try:
        import parmed
        structure = parmed.openmm.load_topology(topology, system, xyz=positions)
        prmtop = os.path.join(out_dir, f"{run_id}_complex.prmtop")
        inpcrd = os.path.join(out_dir, f"{run_id}_complex.inpcrd")
        structure.save(prmtop, overwrite=True)
        structure.save(inpcrd, overwrite=True)
        return prmtop if os.path.exists(prmtop) else None
    except Exception:
        return None


def _smiles_to_sdf(smiles: str, out_path: str) -> bool:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0:
            return False
        AllChem.MMFFOptimizeMolecule(mol)
        Chem.MolToMolFile(mol, out_path)
        return os.path.exists(out_path)
    except Exception:
        return False


def _build_system(protein_pdb: str, ligand_sdf: str):
    """Solvated protein-ligand system; returns (system, topology, positions)."""
    from openmm import unit
    from openmm.app import ForceField, Modeller, PDBFile
    from openff.toolkit import Molecule
    from openmmforcefields.generators import SMIRNOFFTemplateGenerator

    pdb = PDBFile(protein_pdb)
    lig = Molecule.from_file(ligand_sdf)
    ff = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    ff.registerTemplateGenerator(SMIRNOFFTemplateGenerator(molecules=[lig]).generator)
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(ff)
    modeller.addSolvent(ff, padding=1.0 * unit.nanometer,
                        ionicStrength=_SALT_CONC * unit.molar)
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=__import__("openmm").app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=__import__("openmm").app.HBonds,
    )
    return system, modeller.topology, modeller.positions


def _run_simulation(system, topology, positions, out_dir: str, run_id: str) -> dict:
    import openmm
    from openmm import unit
    from openmm.app import DCDReporter, PDBFile, Simulation, StateDataReporter

    system.addForce(openmm.MonteCarloBarostat(
        1.0 * unit.atmosphere, _TEMPERATURE * unit.kelvin))
    integrator = openmm.LangevinMiddleIntegrator(
        _TEMPERATURE * unit.kelvin, 1.0 / unit.picosecond,
        _DT_FS * unit.picoseconds)
    sim = Simulation(topology, system, integrator)
    sim.context.setPositions(positions)
    sim.minimizeEnergy()
    sim.step(_EQUIL_STEPS)

    traj = os.path.join(out_dir, f"{run_id}.dcd")
    log = os.path.join(out_dir, f"{run_id}_energy.csv")
    sim.reporters.append(DCDReporter(traj, _REPORT_INTERVAL))
    sim.reporters.append(StateDataReporter(
        log, _REPORT_INTERVAL,
        step=True, potentialEnergy=True, temperature=True, density=True))
    sim.step(_MD_STEPS)

    final_pdb = os.path.join(out_dir, f"{run_id}_final.pdb")
    state = sim.context.getState(getPositions=True)
    with open(final_pdb, "w") as fh:
        PDBFile.writeFile(topology, state.getPositions(), fh)

    duration_ns = round(_MD_STEPS * _DT_FS / 1000, 1)
    return {"trajectory": traj, "energy_log": log, "final_pdb": final_pdb,
            "duration_ns": duration_ns, "status": "completed"}


def md_simulate_node(state: VTAState) -> VTAState:
    """Run OpenMM MD on each md_candidate; degrade to skip when OpenMM absent."""
    candidates = state.get("md_candidates") or []
    if not candidates:
        state["md_results"] = {}
        return state

    if not _check_openmm():
        state["md_results"] = {c["ligand"]: {"status": "skipped",
                                              "reason": "OpenMM not installed"}
                               for c in candidates}
        state["audit_trail"].append(
            f"[skip] MD-simulate: OpenMM not installed — {len(candidates)} "
            "candidates labelled skipped; install openmm + openmmforcefields")
        state["versions"]["md"] = "skipped (no OpenMM)"
        return state

    out_dir = os.path.join("outputs", state["run_id"], "md")
    os.makedirs(out_dir, exist_ok=True)
    results: dict = {}

    for c in candidates:
        lid, protein = c["ligand"], c["protein"]
        run_id = f"{lid}_{protein}".replace(" ", "_")
        try:
            pdb_path = (state.get("structures") or {}).get(protein, {}).get("pdb_path")
            if not pdb_path:
                raise FileNotFoundError(f"no PDB for {protein}")
            sdf = os.path.join(out_dir, f"{run_id}.sdf")
            if not _smiles_to_sdf(c.get("smiles", ""), sdf):
                raise ValueError(f"SMILES→SDF failed for {lid}")
            sys_, top, pos = _build_system(pdb_path, sdf)
            res = _run_simulation(sys_, top, pos, out_dir, run_id)
            # Emit an Amber topology so md_analyze's MM-GBSA can actually run in
            # production (it reads result["parm"]). Graceful: ParmEd absent or save
            # fails → omit `parm`, MD still completes on RMSD/contacts.
            if _check_parmed():
                parm = _save_amber_topology(sys_, top, pos, out_dir, run_id)
                if parm:
                    res["parm"] = parm
            results[lid] = res
            parm_note = " + prmtop" if res.get("parm") else ""
            state["audit_trail"].append(
                f"MD-simulate[{lid}]: {res['duration_ns']} ns completed → "
                f"{Path(res['trajectory']).name}{parm_note}")
        except Exception as e:
            results[lid] = {"status": "failed", "reason": str(e)}
            state["audit_trail"].append(f"MD-simulate[{lid}]: FAILED — {e}")

    state["md_results"] = results
    dur = _MD_STEPS * _DT_FS / 1000
    state["versions"]["md"] = f"OpenMM ({dur:.0f} ns NPT, AMBER ff14SB + OpenFF 2.0)"
    return state
