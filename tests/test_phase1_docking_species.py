"""Phase-1 docking species, metal, and ensemble seams."""
from __future__ import annotations

from vta.nodes.docking import _dock_mock
from vta.nodes.species_resolution import species_resolution_node
from vta.state import new_state


def _ligand(name="Remdesivir", smiles="CCO"):
    return {
        "name": name,
        "chembl_id": f"CHEMBL_{name}",
        "smiles": smiles,
        "positive_control": True,
    }


def _state():
    st = new_state("x", "phase1")
    st["structures"] = {"PB1": {"pdb_path": "primary.pdb", "ensemble_pdbs": ["alt.pdb"]}}
    st["pockets"] = {"PB1": [{"id": 1, "center": [1, 2, 3], "conservation": 0.8}]}
    return st


def test_species_resolution_node_records_upstream_decisions(monkeypatch):
    import vta.nodes.species_resolution as node

    monkeypatch.setattr(node, "load_ligands", lambda: [_ligand()])
    out = species_resolution_node(new_state("x", "r"))

    assert out["docking_species"]["Remdesivir"]["species_source"] == "parent_surrogate"
    assert "Species resolution: 0 active-form substitutions, 1 parent surrogates" in "\n".join(
        out["audit_trail"])


def test_mock_docking_adds_phase1_species_metal_and_ensemble_fields(monkeypatch):
    import vta.nodes.docking as docking

    monkeypatch.setattr(docking, "load_ligands", lambda: [_ligand()])
    out = _dock_mock(_state())
    row = out["docking_results"][0]

    assert row["species_source"] == "parent_surrogate"
    assert row["dock_species"] == "Remdesivir"
    assert row["metal_model_status"] == "skipped_no_curated_metals"
    assert row["ensemble_size"] == 2
    assert row["ensemble_status"] == "ensemble"
    assert "Docking species: 0 active-form substitutions, 1 parent surrogates" in "\n".join(
        out["audit_trail"])


def test_mock_docking_uses_curated_active_form_when_present(monkeypatch):
    import vta.chem.species as species
    import vta.nodes.docking as docking

    monkeypatch.setitem(species._ACTIVE_SPECIES, "testprodrug", {
        "prodrug": True,
        "active_form": "Test triphosphate",
        "active_form_smiles": "CCN",
        "mechanistic_note": "Test active form.",
    })
    monkeypatch.setattr(docking, "load_ligands", lambda: [_ligand("Testprodrug", "CCO")])
    st = _state()
    st["structures"]["PB1"] = {
        "pdb_path": "primary.pdb",
        "catalytic_metals_by_pocket": {1: [{"element": "Mg", "coord": [0, 0, 0]}]},
    }
    out = _dock_mock(st)
    row = out["docking_results"][0]

    assert row["dock_smiles"] == "CCN"
    assert row["dock_species"] == "Test triphosphate"
    assert row["species_source"] == "curated_active_form"
    assert row["metal_model_status"] == "curated"
