"""chemistry_node — annotation-only chemistry and active-species flags."""
from __future__ import annotations

from vta.chem.filters import chemistry_flags
from vta.chem.species import active_species
from vta.state import VTAState


def chemistry_node(state: VTAState) -> VTAState:
    """Annotate ranked leads with chemistry filters and active species metadata."""
    leads = state.get("lead_candidates") or []
    annotations: dict[str, dict] = {}
    if not leads:
        state["chemistry_annotations"] = annotations
        return state

    prodrugs = flagged = 0
    for lead in leads:
        ligand = lead.get("ligand")
        smiles = lead.get("smiles")
        species = active_species(ligand, smiles)
        flags = chemistry_flags(smiles)
        if species.get("prodrug"):
            prodrugs += 1
        if flags.get("pains") or flags.get("brenk") or flags.get("aggregator"):
            flagged += 1
        lead["active_species"] = species
        lead["chemistry_flags"] = flags
        annotations[ligand] = {"active_species": species, "chemistry_flags": flags}

    state["chemistry_annotations"] = annotations
    state["versions"]["chemistry"] = "RDKit filters + curated active-species map"
    state["audit_trail"].append(
        f"Chemistry: annotated {len(leads)} leads ({prodrugs} prodrug/active-species notes, "
        f"{flagged} chemistry alerts; annotation-only)"
    )
    return state
