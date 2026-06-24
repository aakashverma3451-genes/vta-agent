"""conservation_contacts_node — ligand-contact-weighted conservation (SPEC #4).

SPEC #1 made conservation a per-POCKET scalar, so every ligand docked in the same pocket
gets the SAME value — a uniform offset in `rank.py` that CANNOT reorder leads (the lead
verified post-#1 ranking is byte-identical to the 0.5 baseline). v2 makes it ligand-
specific: score each ligand by the conservation of the residues ITS DOCKED POSE actually
contacts. A ligand gripping the conserved catalytic core (mutation-resistant, durable)
outscores one touching a variable rim — the same per-residue JSD evidence (Capra & Singh
2007), now resolved to the binding interaction.

Mechanism. For each docking record with a saved real pose (`pose_path`+`receptor_path`,
written only by the real Vina path), find the receptor residues with any atom within
R=4 Å of any ligand atom (the same contact definition as md_analyze._contacts_analysis),
look their JSD up in `state["residue_conservation"][protein]` (keyed by resseq, written by
conservation_node), and overwrite `record["conservation"]` with the mean.

Honest fallback (graceful-degradation discipline). No saved pose (mock docking) OR no
residue map (no MSA) → leave the v1 pocket-level value untouched and label the skip. So
mock screens keep the pocket value; real Vina screens get per-ligand values. rank.py is
unchanged — it already reads `record["conservation"]`; v2 just makes that value vary.
"""
from __future__ import annotations

import os

from vta.nodes.conservation import parse_residues
from vta.state import VTAState

R_CONTACT = 4.0     # Å — protein–ligand contact (matches md_analyze._contacts_analysis)


def _parse_ligand_atoms(pdbqt_text: str) -> list[tuple]:
    """(x,y,z) of every atom in the FIRST pose/model of a Vina .pdbqt output."""
    atoms: list[tuple] = []
    for line in pdbqt_text.splitlines():
        if line.startswith("ENDMDL"):
            break                                  # first docked pose only
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            atoms.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        except ValueError:
            continue
    return atoms


def _contact_resseqs(residues: list[dict], ligand_atoms: list[tuple],
                     radius: float = R_CONTACT) -> list[str]:
    """resseq of every receptor residue with an atom within `radius` Å of any ligand atom."""
    r2 = radius * radius
    hits: list[str] = []
    for res in residues:
        for (rx, ry, rz) in res["atoms"]:
            if any((rx - lx) ** 2 + (ry - ly) ** 2 + (rz - lz) ** 2 <= r2
                   for (lx, ly, lz) in ligand_atoms):
                hits.append(res["key"][1])         # resseq, matches residue_conservation
                break
    return hits


def conservation_contacts_node(state: VTAState) -> VTAState:
    rows = state.get("docking_results") or []
    if not rows:
        return state
    res_cons = state.get("residue_conservation") or {}

    receptor_cache: dict[str, list[dict]] = {}     # parse each receptor PDBQT once
    weighted = no_pose = no_map = 0
    for r in rows:
        pose, receptor = r.get("pose_path"), r.get("receptor_path")
        if not (pose and receptor and os.path.exists(pose) and os.path.exists(receptor)):
            no_pose += 1
            continue                               # mock docking → keep v1 pocket value
        cons_map = res_cons.get(r.get("protein"))
        if not cons_map:
            no_map += 1
            continue                               # no MSA → keep v1 pocket value
        residues = receptor_cache.get(receptor)
        if residues is None:
            try:
                with open(receptor) as fh:
                    residues = parse_residues(fh.read())
            except OSError:
                residues = []
            receptor_cache[receptor] = residues
        try:
            with open(pose) as fh:
                lig_atoms = _parse_ligand_atoms(fh.read())
        except OSError:
            lig_atoms = []
        if not (residues and lig_atoms):
            no_pose += 1
            continue
        contacts = _contact_resseqs(residues, lig_atoms)
        scores = [cons_map[k] for k in contacts if k in cons_map]
        if scores:
            r["conservation"] = round(sum(scores) / len(scores), 3)
            weighted += 1

    if weighted:
        state["audit_trail"].append(
            f"contact-conservation: ligand-weighted {weighted} poses "
            f"(per-residue JSD over {R_CONTACT}Å contacts; ranking now ligand-specific)")
        state["versions"]["conservation_contacts"] = "ligand-contact-weighted JSD (4Å)"
    if no_pose:
        state["audit_trail"].append(
            f"[skip] contact-conservation: no pose for {no_pose} records "
            "(kept pocket-level conservation)")
    if no_map:
        state["audit_trail"].append(
            f"[skip] contact-conservation: no residue map for {no_map} records "
            "(no MSA; kept pocket-level conservation)")
    return state
