"""Chemistry annotation tests — hermetic."""
from __future__ import annotations

from vta.chem.filters import chemistry_flags
from vta.chem.species import active_species, resolve_docking_species
from vta.nodes.chemistry import chemistry_node
from vta.provenance import input_hash, score_provenance
from vta.state import new_state


def test_active_species_maps_committed_controls():
    rem = active_species("Remdesivir", "parent")
    assert rem["prodrug"] is True
    assert "triphosphate" in rem["active_form"].lower()
    mol = active_species("Molnupiravir")
    assert mol["prodrug"] is True
    assert "triphosphate" in mol["active_form"].lower()
    rib = active_species("Ribavirin")
    assert rib["prodrug"] is False
    assert "triphosphate" in rib["active_form"].lower()


def test_active_species_unknown_uses_parent():
    out = active_species("Lopinavir", "CC")
    assert out["prodrug"] is False
    assert out["active_form"] == "Lopinavir"
    assert out["active_form_smiles"] == "CC"


def test_resolve_docking_species_labels_parent_surrogate():
    out = resolve_docking_species({"name": "Remdesivir", "smiles": "parent"})
    assert out["dock_smiles"] == "parent"
    assert out["species_source"] == "parent_surrogate"
    assert out["uses_active_form"] is False
    assert "no curated active-form SMILES" in out["mechanistic_caveat"]


def test_resolve_docking_species_uses_curated_active_form(monkeypatch):
    import vta.chem.species as species

    monkeypatch.setitem(species._ACTIVE_SPECIES, "testprodrug", {
        "prodrug": True,
        "active_form": "Test triphosphate",
        "active_form_smiles": "CCN",
        "mechanistic_note": "Test active form.",
    })
    out = resolve_docking_species({"name": "Testprodrug", "smiles": "CCO"})
    assert out["dock_smiles"] == "CCN"
    assert out["dock_species"] == "Test triphosphate"
    assert out["species_source"] == "curated_active_form"
    assert out["uses_active_form"] is True


def test_chemistry_flags_clean_molecule():
    flags = chemistry_flags("CCO")
    assert flags["pains"] == []
    assert flags["aggregator"] is False
    assert flags["beyond_ro5"] is False


def test_chemistry_flags_known_alert():
    flags = chemistry_flags("O=C1NC(=S)SC1=Cc1ccc(O)cc1")
    assert flags["pains"] or flags["brenk"]


def test_chemistry_node_annotates_leads():
    st = new_state("x", "x")
    st["lead_candidates"] = [
        {"ligand": "Remdesivir", "smiles": "CCO", "score": 0.9},
        {"ligand": "Lopinavir", "smiles": "CC", "score": 0.8},
    ]
    out = chemistry_node(st)
    assert out["lead_candidates"][0]["active_species"]["prodrug"] is True
    assert out["lead_candidates"][1]["active_species"]["active_form"] == "Lopinavir"
    assert "Chemistry: annotated 2 leads" in "\n".join(out["audit_trail"])


def test_score_provenance_is_stable():
    a = input_hash({"b": 2, "a": 1})
    b = input_hash({"a": 1, "b": 2})
    assert a == b
    prov = score_provenance("tool", "1.0", seed=42, inputs={"x": 1})
    assert prov["tool"] == "tool"
    assert prov["version"] == "1.0"
    assert prov["seed"] == 42
    assert len(prov["input_hash"]) == 16
