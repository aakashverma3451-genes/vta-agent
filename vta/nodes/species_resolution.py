"""species_resolution_node — choose the ligand species used by docking."""
from __future__ import annotations

from vta.chem.species import resolve_docking_species
from vta.data.ligands import load_ligands
from vta.state import VTAState


def species_resolution_node(state: VTAState) -> VTAState:
    """Resolve active-form/prodrug handling before the docking node."""
    ligands = load_ligands()
    resolved = {lig["name"]: resolve_docking_species(lig) for lig in ligands}
    active = sum(1 for r in resolved.values() if r.get("uses_active_form"))
    surrogates = sum(1 for r in resolved.values() if r.get("species_source") == "parent_surrogate")

    state["docking_species"] = resolved
    state["versions"]["docking_species"] = "curated active-species resolver"
    state["audit_trail"].append(
        f"Species resolution: {active} active-form substitutions, "
        f"{surrogates} parent surrogates before docking"
    )
    return state
