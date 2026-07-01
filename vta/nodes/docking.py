"""docking_node — virtual screening with real Vina or deterministic fallback.

The node resolves active/prodrug docking species before ligand prep, uses real Vina
when Vina plus RDKit/Meeko are available, and otherwise emits same-shaped mock rows.
Ranking downstream is unchanged.
"""
from __future__ import annotations

import hashlib
import os
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path

from vta.data.ligands import load_ligands
from vta.docking.seams import ensemble_conformers, ligand_for_docking, metal_model, species_fields
from vta.provenance import score_provenance
from vta.state import VTAState

_TOP_POCKETS_MOCK = 3
_REAL_TOP_POCKETS = 1          # real Vina is expensive (emulated) → dock the top pocket
_BOX = 22.0                    # docking box edge, Å
_EXHAUSTIVENESS = 8
_LIG_CACHE = "structures/.ligand_pdbqt"   # prepped ligands reused across runs


def _seed(*parts: str) -> int:
    return int(hashlib.md5("|".join(parts).encode()).hexdigest()[:8], 16)


# ── real Vina path ───────────────────────────────────────────────────────────
def _vina_bin() -> str | None:
    from vta.toolconfig import find_tool
    return find_tool("vina", "VINA_BIN")


def _prep_available() -> bool:
    try:
        import meeko  # noqa: F401
        import rdkit  # noqa: F401
        return True
    except Exception:
        return False


def _heavy_atoms(smiles: str | None, default: int = 30) -> int:
    """Real heavy-atom count from SMILES (for honest ligand efficiency)."""
    if not smiles:
        return default
    try:
        from rdkit import Chem
        m = Chem.MolFromSmiles(smiles)
        return m.GetNumHeavyAtoms() if m else default
    except Exception:
        return default


def _prep_ligand(lig: dict, cache_dir: str) -> str | None:
    """SMILES → 3D (RDKit ETKDG) → PDBQT (Meeko). Cached by ChEMBL id."""
    if not lig.get("smiles"):
        return None
    os.makedirs(cache_dir, exist_ok=True)
    out = os.path.join(cache_dir, f"{lig.get('dock_cache_id') or lig['chembl_id']}.pdbqt")
    if os.path.exists(out):
        return out
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.MolFromSmiles(lig["smiles"])
    if m is None:
        return None
    m = Chem.AddHs(m)
    if AllChem.EmbedMolecule(m, AllChem.ETKDGv3()) != 0:
        return None
    AllChem.MMFFOptimizeMolecule(m)
    from meeko import MoleculePreparation
    try:                                   # non-deprecated path (Meeko ≥0.5)
        from meeko import PDBQTWriterLegacy
        setups = MoleculePreparation().prepare(m)
        pdbqt, ok, _ = PDBQTWriterLegacy.write_string(setups[0])
        if not ok:
            return None
    except Exception:                      # older Meeko fallback
        prep = MoleculePreparation(); prep.prepare(m)
        pdbqt = prep.write_pdbqt_string()
    with open(out, "w") as fh:
        fh.write(pdbqt)
    return out


def _prep_receptor(pdb_path: str, out_prefix: str) -> str | None:
    """Protein PDB → PDBQT via Meeko's mk_prepare_receptor CLI.

    Uses sys.executable + the tool that ships beside it, so this works without the
    venv being on PATH — required for the autonomous `vta run` (no PATH setup).
    """
    pdbqt = f"{out_prefix}.pdbqt"
    if os.path.exists(pdbqt):
        return pdbqt
    bindir = Path(sys.executable).parent
    tool = (shutil.which("mk_prepare_receptor.py")
            or (str(bindir / "mk_prepare_receptor.py")
                if (bindir / "mk_prepare_receptor.py").exists()
                else "mk_prepare_receptor.py"))
    subprocess.run(
        [sys.executable, tool, "--read_pdb", pdb_path, "-o", out_prefix, "-p",
         "--allow_bad_res", "--default_altloc", "A"],
        capture_output=True, timeout=600,
    )
    if os.path.exists(pdbqt):
        return pdbqt
    # Meeko-first failed (e.g. the deterministic Meeko 0.7.1 Mpro residue-template H bug).
    # Fall back to OpenBabel, which prepares a rigid receptor PDBQT from many structures
    # Meeko refuses. Meeko remains primary, so structures it handles (e.g. TiLV 8PSO) are
    # unaffected; the fallback only runs when Meeko produced nothing.
    return _prep_receptor_obabel(pdb_path, pdbqt)


def _prep_receptor_obabel(pdb_path: str, pdbqt: str) -> str | None:
    """Fallback receptor prep via OpenBabel → rigid PDBQT (Gasteiger charges, AD4 types).

    Strips HETATM (waters / co-crystal ligand / ions) so the rigid receptor excludes any
    bound ligand, adds hydrogens, and writes a rigid (-xr) receptor PDBQT. A REMARK records
    that the fallback was used so downstream provenance can report the prep method honestly.
    """
    try:
        from openbabel import pybel
    except Exception:
        return None
    protein = [ln for ln in open(pdb_path) if ln.startswith(("ATOM", "TER"))]
    if not protein:
        return None
    tmp_pdb = f"{pdbqt}.receptor.pdb"
    with open(tmp_pdb, "w") as fh:
        fh.write("".join(protein))
        fh.write("END\n")
    try:
        mol = next(pybel.readfile("pdb", tmp_pdb))
    except Exception:
        return None
    try:
        mol.addh()                       # add hydrogens for protonation
        mol.write("pdbqt", pdbqt, opt={"r": True}, overwrite=True)  # -xr: rigid receptor
    except Exception:
        return None
    if not os.path.exists(pdbqt) or os.path.getsize(pdbqt) == 0:
        return None
    # Tag provenance: mark the receptor as OpenBabel-prepared.
    body = open(pdbqt).read()
    with open(pdbqt, "w") as fh:
        fh.write("REMARK VTA receptor prep: OpenBabel fallback (Meeko declined)\n")
        fh.write(body)
    return pdbqt


def _run_vina(vina: str, receptor: str, ligand: str, center: list, out: str) -> float | None:
    """Run Vina and return the best-mode affinity (kcal/mol), or None on failure."""
    p = subprocess.run(
        [vina, "--receptor", receptor, "--ligand", ligand,
         "--center_x", str(center[0]), "--center_y", str(center[1]),
         "--center_z", str(center[2]),
         "--size_x", str(_BOX), "--size_y", str(_BOX), "--size_z", str(_BOX),
         "--exhaustiveness", str(_EXHAUSTIVENESS), "--num_modes", "5",
         "--seed", "42", "--out", out],   # fixed seed → reproducible ΔG
        capture_output=True, text=True, timeout=900,
    )
    for line in p.stdout.splitlines():
        m = re.match(r"\s*1\s+(-?\d+\.\d+)", line)      # mode 1 = best
        if m:
            return float(m.group(1))
    return None


def _dock_real(state: VTAState, vina: str) -> VTAState:
    ligands = load_ligands()
    results, n_docked = [], 0
    active_substitutions = parent_surrogates = metal_skips = ensemble_jobs = 0
    for protein, plist in (state.get("pockets") or {}).items():
        struct = (state.get("structures") or {}).get(protein, {})
        if not struct.get("pdb_path") or not plist:
            continue
        conformers = ensemble_conformers(struct)
        receptors = []
        for conformer in conformers:
            prefix = os.path.join("structures", f"{state['run_id']}_{protein}_{conformer['label']}_rec")
            receptor = _prep_receptor(conformer["pdb_path"], prefix)
            if receptor:
                receptors.append({**conformer, "receptor": receptor})
        if not receptors:
            state["audit_trail"].append(f"Vina[{protein}]: receptor prep failed; skipped")
            continue
        if len(receptors) > 1:
            ensemble_jobs += 1
        for pocket in plist[:_REAL_TOP_POCKETS]:
            metal = metal_model(struct, pocket)
            if metal["status"] != "curated":
                metal_skips += 1
            for lig in ligands:
                # Per-ligand isolation: one bad SMILES / failed dock skips that ligand,
                # it does NOT abort the whole real screen (which would fall back to mock).
                try:
                    dock_lig, resolution = ligand_for_docking(lig)
                    if resolution["uses_active_form"]:
                        active_substitutions += 1
                    elif resolution["species_source"] == "parent_surrogate":
                        parent_surrogates += 1
                    lig_pdbqt = _prep_ligand(dock_lig, _LIG_CACHE)
                    if not lig_pdbqt:
                        state["audit_trail"].append(
                            f"Vina[{protein}]: {lig['name']} ligand prep failed; skipped")
                        continue
                    best = None
                    for rec in receptors:
                        out = os.path.join(
                            "structures",
                            f"{state['run_id']}_{protein}_{lig['chembl_id']}_p{pocket['id']}_{rec['label']}.pdbqt")
                        dG = _run_vina(vina, rec["receptor"], lig_pdbqt, pocket["center"], out)
                        if dG is not None and (best is None or dG < best["dG"]):
                            best = {"dG": dG, "pose": out, "receptor": rec["receptor"],
                                    "conformer": rec["label"]}
                    dG = best["dG"] if best else None
                    if dG is None:
                        state["audit_trail"].append(
                            f"Vina[{protein}]: {lig['name']} produced no pose; skipped")
                        continue
                    heavy = _heavy_atoms(dock_lig.get("smiles"))
                    results.append({
                        "protein": protein, "pocket": pocket["id"],
                        "ligand": lig["name"], "ligand_id": lig["chembl_id"],
                        "smiles": lig["smiles"], "positive_control": lig["positive_control"],
                        "dG": dG, "rmsd": 0.0, "le": round(dG / heavy, 3),
                        "heavy_atoms": heavy, "conservation": pocket["conservation"],
                        **species_fields(resolution, metal, len(receptors)),
                        "ensemble_conformer": best["conformer"],
                        "dG_provenance": score_provenance(
                            "AutoDock Vina", "1.2.5", seed=42, inputs={
                                "protein": protein, "pocket": pocket["id"],
                                "ligand_id": lig["chembl_id"],
                                "dock_species": resolution.get("dock_species"),
                                "species_source": resolution.get("species_source"),
                                "center": pocket["center"],
                                "box": _BOX, "exhaustiveness": _EXHAUSTIVENESS,
                            }),
                        # saved pose + receptor let the DL-rescore seam (rescore_node)
                        # re-score this pose later without re-docking.
                        "pose_path": best["pose"], "receptor_path": best["receptor"],
                    })
                    n_docked += 1
                except Exception as e:
                    state["audit_trail"].append(
                        f"Vina[{protein}]: {lig['name']} error ({type(e).__name__}: {e}); skipped")
    state["docking_results"] = results
    state["audit_trail"].append(f"Vina: {n_docked} real docks over {len(ligands)} ligands")
    state["audit_trail"].append(
        f"Docking species: {active_substitutions} active-form substitutions, "
        f"{parent_surrogates} parent surrogates")
    if metal_skips:
        state["audit_trail"].append(
            f"[skip] Metals: {metal_skips} pocket(s) lacked curated catalytic Mg/Mn coordinates")
    if ensemble_jobs == 0:
        state["audit_trail"].append("[skip] Ensemble docking: no receptor ensemble supplied")
    else:
        state["audit_trail"].append(f"Ensemble docking: used ensembles for {ensemble_jobs} protein(s)")
    state["versions"]["docking"] = "AutoDock Vina 1.2.5"
    state["versions"]["ligands"] = f"ChEMBL x{len(ligands)}"
    state["versions"]["docking_species"] = "curated active-species resolver"
    return state

# ── mock fallback (real ligands, fabricated scores) ──────────────────────────
def _dock_mock(state: VTAState) -> VTAState:
    ligands = load_ligands()
    results = []
    active_substitutions = parent_surrogates = metal_skips = ensemble_jobs = 0
    for protein, plist in (state.get("pockets") or {}).items():
        struct = (state.get("structures") or {}).get(protein, {})
        conformers = ensemble_conformers(struct)
        ensemble_size = max(1, len(conformers))
        if len(conformers) > 1:
            ensemble_jobs += 1
        for pocket in plist[:_TOP_POCKETS_MOCK]:
            metal = metal_model(struct, pocket)
            if metal["status"] != "curated":
                metal_skips += 1
            for lig in ligands:
                dock_lig, resolution = ligand_for_docking(lig)
                if resolution["uses_active_form"]:
                    active_substitutions += 1
                elif resolution["species_source"] == "parent_surrogate":
                    parent_surrogates += 1
                rng = random.Random(
                    _seed(state["run_id"], protein, str(pocket["id"]), lig["chembl_id"]))
                dG = round(rng.uniform(-9.0, -5.0), 2)
                rmsd = round(rng.uniform(0.5, 3.0), 2)
                heavy = _heavy_atoms(dock_lig.get("smiles"))
                results.append({
                    "protein": protein, "pocket": pocket["id"],
                    "ligand": lig["name"], "ligand_id": lig["chembl_id"],
                    "smiles": lig["smiles"], "positive_control": lig["positive_control"],
                    "dG": dG, "rmsd": rmsd, "le": round(dG / heavy, 3),
                    "heavy_atoms": heavy, "conservation": pocket["conservation"],
                    **species_fields(resolution, metal, ensemble_size),
                    "dG_provenance": score_provenance(
                        "MOCK-vina-shaped", "deterministic", seed=_seed(
                            state["run_id"], protein, str(pocket["id"]), lig["chembl_id"]),
                        inputs={"protein": protein, "pocket": pocket["id"],
                                "ligand_id": lig["chembl_id"],
                                "dock_species": resolution.get("dock_species"),
                                "species_source": resolution.get("species_source")}),
                })
    state["docking_results"] = results
    state["audit_trail"].append(
        f"[MOCK] Docking: {len(results)} pairs over {len(ligands)} real ligands")
    state["audit_trail"].append(
        f"Docking species: {active_substitutions} active-form substitutions, "
        f"{parent_surrogates} parent surrogates")
    if metal_skips:
        state["audit_trail"].append(
            f"[skip] Metals: {metal_skips} pocket(s) lacked curated catalytic Mg/Mn coordinates")
    if ensemble_jobs == 0:
        state["audit_trail"].append("[skip] Ensemble docking: no receptor ensemble supplied")
    else:
        state["audit_trail"].append(f"Ensemble docking: detected ensembles for {ensemble_jobs} protein(s)")
    state["versions"]["docking"] = "MOCK-vina-shaped over ChEMBL ligands"
    state["versions"]["ligands"] = f"ChEMBL x{len(ligands)}"
    state["versions"]["docking_species"] = "curated active-species resolver"
    return state

# ── node ─────────────────────────────────────────────────────────────────────
def docking_node(state: VTAState) -> VTAState:
    """Real Vina if engine + prep stack available, else mock."""
    vina = _vina_bin()
    if vina and _prep_available():
        try:
            return _dock_real(state, vina)
        except Exception as e:
            state["audit_trail"].append(f"Vina: WARNING real docking failed ({e}); mock")
    return _dock_mock(state)


# back-compat alias used by the graph / tests
docking_node_mock = docking_node
