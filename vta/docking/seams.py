"""Phase-1 docking seams for species, metals, and receptor ensembles."""
from __future__ import annotations

import hashlib
from typing import Any

from vta.chem.species import resolve_docking_species


def ligand_for_docking(lig: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Return a ligand copy with the resolved docking species SMILES."""
    resolution = resolve_docking_species(lig)
    out = dict(lig)
    out["smiles"] = resolution.get("dock_smiles")
    if resolution.get("uses_active_form"):
        sid = hashlib.md5((resolution.get("dock_smiles") or "").encode()).hexdigest()[:10]
        out["dock_cache_id"] = f"{lig['chembl_id']}_{sid}"
    else:
        out["dock_cache_id"] = lig["chembl_id"]
    return out, resolution


def metal_model(struct: dict[str, Any], pocket: dict[str, Any]) -> dict[str, Any]:
    """Return a curated catalytic-metal model, or a labelled skip."""
    pocket_id = str(pocket.get("id"))
    by_pocket = struct.get("catalytic_metals_by_pocket") or {}
    metals = (
        pocket.get("catalytic_metals")
        or pocket.get("metals")
        or by_pocket.get(pocket_id)
        or by_pocket.get(pocket.get("id"))
    )
    if metals:
        return {"status": "curated", "metals": metals}
    return {"status": "skipped_no_curated_metals", "metals": []}


def ensemble_conformers(struct: dict[str, Any]) -> list[dict[str, Any]]:
    """Return primary plus any explicitly supplied receptor conformers."""
    primary = struct.get("pdb_path")
    if not primary:
        return []
    conformers = [{"label": "primary", "pdb_path": primary}]
    for i, path in enumerate(struct.get("ensemble_pdbs") or struct.get("receptor_ensemble") or [], 1):
        if path and path != primary:
            conformers.append({"label": f"ensemble_{i}", "pdb_path": path})
    return conformers


def species_fields(resolution: dict[str, Any], metal: dict[str, Any], ensemble_size: int) -> dict[str, Any]:
    return {
        "dock_smiles": resolution.get("dock_smiles"),
        "dock_species": resolution.get("dock_species"),
        "species_source": resolution.get("species_source"),
        "mechanistic_caveat": resolution.get("mechanistic_caveat"),
        "active_species": resolution.get("active_species"),
        "metal_model_status": metal.get("status"),
        "metal_model": metal,
        "ensemble_size": ensemble_size,
        "ensemble_status": "ensemble" if ensemble_size > 1 else "single_receptor",
    }
