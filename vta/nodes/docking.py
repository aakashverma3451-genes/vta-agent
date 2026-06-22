"""docking_node — virtual screening.

Phase 2c: REAL AutoDock Vina when available, else a real-shaped MOCK.

`docking_node` auto-detects the Vina engine (`VINA_BIN` or `vina` on PATH) plus the
RDKit/Meeko prep stack. If all are present it docks for real:

    SMILES --RDKit 3D--> --Meeko--> ligand.pdbqt
    receptor.pdb --Meeko mk_prepare_receptor--> receptor.pdbqt
    vina --center (from FPocket pocket) --> best affinity (ΔG, kcal/mol)

…otherwise it degrades to the deterministic mock (same don't-crash discipline as the
other real nodes). Either way the records have the same shape, so ranking downstream
is unchanged. Real Vina under x86 emulation is slow, so the real path docks the top
pocket per protein (bounded), with prepped ligand PDBQTs cached across runs.
"""
from __future__ import annotations

import hashlib
import os
import random
import re
import shutil
import subprocess

from vta.data.ligands import load_ligands
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
    return os.environ.get("VINA_BIN") or shutil.which("vina")


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
    out = os.path.join(cache_dir, f"{lig['chembl_id']}.pdbqt")
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
    """Protein PDB → PDBQT via Meeko's mk_prepare_receptor CLI."""
    pdbqt = f"{out_prefix}.pdbqt"
    if os.path.exists(pdbqt):
        return pdbqt
    tool = shutil.which("mk_prepare_receptor.py") or "mk_prepare_receptor.py"
    r = subprocess.run(
        ["python", tool, "--read_pdb", pdb_path, "-o", out_prefix, "-p",
         "--allow_bad_res", "--default_altloc", "A"],
        capture_output=True, timeout=600,
    )
    return pdbqt if os.path.exists(pdbqt) else None


def _run_vina(vina: str, receptor: str, ligand: str, center: list, out: str) -> float | None:
    """Run Vina and return the best-mode affinity (kcal/mol), or None on failure."""
    p = subprocess.run(
        [vina, "--receptor", receptor, "--ligand", ligand,
         "--center_x", str(center[0]), "--center_y", str(center[1]),
         "--center_z", str(center[2]),
         "--size_x", str(_BOX), "--size_y", str(_BOX), "--size_z", str(_BOX),
         "--exhaustiveness", str(_EXHAUSTIVENESS), "--num_modes", "5", "--out", out],
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
    for protein, plist in (state.get("pockets") or {}).items():
        struct = (state.get("structures") or {}).get(protein, {})
        if not struct.get("pdb_path") or not plist:
            continue
        receptor = _prep_receptor(struct["pdb_path"],
                                  os.path.join("structures", f"{state['run_id']}_{protein}_rec"))
        if not receptor:
            state["audit_trail"].append(f"Vina[{protein}]: receptor prep failed; skipped")
            continue
        for pocket in plist[:_REAL_TOP_POCKETS]:
            for lig in ligands:
                # Per-ligand isolation: one bad SMILES / failed dock skips that ligand,
                # it does NOT abort the whole real screen (which would fall back to mock).
                try:
                    lig_pdbqt = _prep_ligand(lig, _LIG_CACHE)
                    if not lig_pdbqt:
                        state["audit_trail"].append(
                            f"Vina[{protein}]: {lig['name']} ligand prep failed; skipped")
                        continue
                    out = os.path.join("structures",
                                       f"{state['run_id']}_{protein}_{lig['chembl_id']}_p{pocket['id']}.pdbqt")
                    dG = _run_vina(vina, receptor, lig_pdbqt, pocket["center"], out)
                    if dG is None:
                        state["audit_trail"].append(
                            f"Vina[{protein}]: {lig['name']} produced no pose; skipped")
                        continue
                    heavy = _heavy_atoms(lig.get("smiles"))
                    results.append({
                        "protein": protein, "pocket": pocket["id"],
                        "ligand": lig["name"], "ligand_id": lig["chembl_id"],
                        "smiles": lig["smiles"], "positive_control": lig["positive_control"],
                        "dG": dG, "rmsd": 0.0, "le": round(dG / heavy, 3),
                        "heavy_atoms": heavy, "conservation": pocket["conservation"],
                    })
                    n_docked += 1
                except Exception as e:
                    state["audit_trail"].append(
                        f"Vina[{protein}]: {lig['name']} error ({type(e).__name__}: {e}); skipped")
    state["docking_results"] = results
    state["audit_trail"].append(f"Vina: {n_docked} real docks over {len(ligands)} ligands")
    state["versions"]["docking"] = "AutoDock Vina 1.2.5"
    state["versions"]["ligands"] = f"ChEMBL x{len(ligands)}"
    return state


# ── mock fallback (real ligands, fabricated scores) ──────────────────────────
def _dock_mock(state: VTAState) -> VTAState:
    ligands = load_ligands()
    results = []
    for protein, plist in (state.get("pockets") or {}).items():
        for pocket in plist[:_TOP_POCKETS_MOCK]:
            for lig in ligands:
                rng = random.Random(
                    _seed(state["run_id"], protein, str(pocket["id"]), lig["chembl_id"]))
                dG = round(rng.uniform(-9.0, -5.0), 2)
                rmsd = round(rng.uniform(0.5, 3.0), 2)
                heavy = _heavy_atoms(lig.get("smiles"))
                results.append({
                    "protein": protein, "pocket": pocket["id"],
                    "ligand": lig["name"], "ligand_id": lig["chembl_id"],
                    "smiles": lig["smiles"], "positive_control": lig["positive_control"],
                    "dG": dG, "rmsd": rmsd, "le": round(dG / heavy, 3),
                    "heavy_atoms": heavy, "conservation": pocket["conservation"],
                })
    state["docking_results"] = results
    state["audit_trail"].append(
        f"[MOCK] Docking: {len(results)} pairs over {len(ligands)} real ligands")
    state["versions"]["docking"] = "MOCK-vina-shaped over ChEMBL ligands"
    state["versions"]["ligands"] = f"ChEMBL x{len(ligands)}"
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
