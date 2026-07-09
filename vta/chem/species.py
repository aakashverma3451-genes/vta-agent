"""Active-species and prodrug annotations for antiviral ligands."""
from __future__ import annotations

from typing import Any


_ACTIVE_SPECIES: dict[str, dict[str, Any]] = {
    "remdesivir": {
        "prodrug": True,
        "active_form": "Remdesivir triphosphate (GS-443902)",
        "active_form_smiles": None,
        "mechanistic_note": "Parent prodrug; RdRp active species is intracellular triphosphate.",
    },
    "sofosbuvir": {
        "prodrug": True,
        "active_form": "Sofosbuvir triphosphate (GS-461203)",
        "active_form_smiles": None,
        "mechanistic_note": "Parent prodrug; NS5B/RdRp active species is triphosphate.",
    },
    "molnupiravir": {
        "prodrug": True,
        "active_form": "NHC triphosphate",
        "active_form_smiles": None,
        "mechanistic_note": "Prodrug of beta-D-N4-hydroxycytidine; active species is triphosphate.",
    },
    "ribavirin": {
        "prodrug": False,
        "active_form": "Ribavirin triphosphate",
        "active_form_smiles": None,
        "mechanistic_note": "Nucleoside analog; triphosphate is the polymerase-relevant species.",
    },
    "favipiravir": {
        "prodrug": False,
        "active_form": "Favipiravir ribofuranosyl-5'-triphosphate",
        "active_form_smiles": None,
        "mechanistic_note": "Base analog; ribosylated triphosphate is polymerase-relevant.",
    },
    "tenofovir": {
        "prodrug": False,
        "active_form": "Tenofovir diphosphate",
        "active_form_smiles": None,
        "mechanistic_note": "Nucleotide analog; diphosphate is the active chain terminator.",
    },
    "lamivudine": {
        "prodrug": False,
        "active_form": "Lamivudine triphosphate",
        "active_form_smiles": None,
        "mechanistic_note": "Nucleoside analog; triphosphate is the active species.",
    },
    "entecavir": {
        "prodrug": False,
        "active_form": "Entecavir triphosphate",
        "active_form_smiles": None,
        "mechanistic_note": "Guanosine analog; triphosphate is the active species.",
    },
    "acyclovir": {
        "prodrug": False,
        "active_form": "Acyclovir triphosphate",
        "active_form_smiles": None,
        "mechanistic_note": "Nucleoside analog; triphosphate is the active species.",
    },
    "ganciclovir": {
        "prodrug": False,
        "active_form": "Ganciclovir triphosphate",
        "active_form_smiles": None,
        "mechanistic_note": "Nucleoside analog; triphosphate is the active species.",
    },
    "valacyclovir": {
        "prodrug": True,
        "active_form": "Acyclovir triphosphate",
        "active_form_smiles": None,
        "mechanistic_note": "Prodrug of acyclovir; triphosphate is the active species.",
    },
    "zidovudine": {
        "prodrug": False,
        "active_form": "Zidovudine triphosphate",
        "active_form_smiles": None,
        "mechanistic_note": "Nucleoside analog; triphosphate is the active species.",
    },
    "emtricitabine": {
        "prodrug": False,
        "active_form": "Emtricitabine triphosphate",
        "active_form_smiles": None,
        "mechanistic_note": "Nucleoside analog; triphosphate is the active species.",
    },
    "abacavir": {
        "prodrug": True,
        "active_form": "Carbovir triphosphate",
        "active_form_smiles": None,
        "mechanistic_note": "Converted intracellularly to carbovir triphosphate.",
    },
}


def active_species(ligand_name: str | None, smiles: str | None = None) -> dict[str, Any]:
    """Return active-form annotation for a ligand name."""
    key = (ligand_name or "").strip().lower()
    rec = _ACTIVE_SPECIES.get(key)
    if rec:
        return {"parent": ligand_name, "parent_smiles": smiles, **rec}
    return {
        "parent": ligand_name,
        "parent_smiles": smiles,
        "prodrug": False,
        "active_form": ligand_name,
        "active_form_smiles": smiles,
        "mechanistic_note": "No curated active-species transform; parent used as screening species.",
    }


def resolve_docking_species(ligand: dict[str, Any]) -> dict[str, Any]:
    """Resolve the chemical species that should be prepared for docking.

    The active-form table is authoritative only when it contains a curated active-form
    SMILES. If a ligand is known to require an active metabolite but the structure is
    not curated yet, docking stays on the parent and the record is labelled as a
    parent surrogate. This avoids silently fabricating triphosphate structures.
    """
    name = ligand.get("name") or ligand.get("ligand")
    parent_smiles = ligand.get("smiles")
    species = active_species(name, parent_smiles)
    active_smiles = species.get("active_form_smiles")
    uses_active = bool(active_smiles and active_smiles != parent_smiles)

    if uses_active:
        source = "curated_active_form"
        dock_smiles = active_smiles
        dock_species = species.get("active_form") or name
        caveat = ""
    elif species.get("active_form") != name and not active_smiles:
        source = "parent_surrogate"
        dock_smiles = parent_smiles
        dock_species = name
        caveat = (
            "Active form is known but no curated active-form SMILES is committed; "
            "parent structure was docked as a labelled surrogate."
        )
    else:
        source = "parent"
        dock_smiles = parent_smiles
        dock_species = name
        caveat = ""

    return {
        "parent": name,
        "parent_smiles": parent_smiles,
        "dock_species": dock_species,
        "dock_smiles": dock_smiles,
        "species_source": source,
        "uses_active_form": uses_active,
        "mechanistic_caveat": caveat,
        "active_species": species,
    }
