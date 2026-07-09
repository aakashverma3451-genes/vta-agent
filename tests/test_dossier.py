"""Tests for the DossierBuilder node (Phase R / R1). Hermetic — no network, no docking."""
from __future__ import annotations

from vta.nodes.dossier import dossier_node
from vta.state import new_state


def _state_with(structures: dict, pockets: dict | None = None) -> dict:
    st = new_state("dossier-t", "dossier-t")
    st["structures"] = structures
    if pockets is not None:
        st["pockets"] = pockets
    return st


def test_dossier_grades_the_three_wired_targets():
    st = _state_with(
        {"MPRO": {"method": "experimental", "source": "7L11:A", "mean_plddt": None},
         "NS5B": {"method": "experimental", "source": "2XI3:A", "mean_plddt": None},
         "PB1": {"method": "experimental", "source": "8PSO:B", "mean_plddt": None}},
    )
    d = dossier_node(st)["target_dossier"]
    # Mpro: protease, non-metal, powered
    assert d["MPRO"]["target_class"] == "protease"
    assert d["MPRO"]["metal_dependence"] is False
    assert d["MPRO"]["benchmarkability"]["status"] == "powered"
    # NS5B (HCV): polymerase, metal-dependent biology, un_benchmarkable per frozen benchmark
    assert d["NS5B"]["metal_dependence"] is True
    assert d["NS5B"]["benchmarkability"]["status"] == "un_benchmarkable"
    # PB1 (TiLV): metal-dependent polymerase biology, but only UNDERPOWERED (still dockable)
    assert d["PB1"]["metal_dependence"] is True
    assert d["PB1"]["benchmarkability"]["status"] == "underpowered"


def test_dossier_records_provenance_and_pocket_metal():
    st = _state_with(
        {"MPRO": {"method": "experimental", "source": "7L11:A", "mean_plddt": None}},
        {"MPRO": [{"id": 1, "center": [1, 2, 3], "volume": 700.0, "druggability": 0.8,
                   "metal_ions": ["MG"]}]},
    )
    d = dossier_node(st)["target_dossier"]["MPRO"]
    assert d["structure_quality"]["provenance"] == "experimental"
    assert d["pocket_descriptors"]["available"] is True
    assert d["pocket_descriptors"]["metal_ions_present"] is True   # metal in THIS pocket


def test_dossier_flags_predicted_low_plddt_and_unknown_benchmark():
    st = _state_with(
        {"GP120": {"method": "esmfold", "mean_plddt": 55.0,
                   "structure_qc": {"status": "flagged", "issues": ["low_plddt"]}}},
    )
    d = dossier_node(st)["target_dossier"]["GP120"]
    assert d["structure_quality"]["provenance"] == "predicted"
    assert d["structure_quality"]["qc_status"] == "flagged"
    assert d["structure_quality"]["binding_site_plddt"] == 55.0
    assert d["benchmarkability"]["status"] == "unknown"            # no frozen benchmark covers it


def test_dossier_never_fabricates_missing_fields():
    st = _state_with({"ORF9": {"method": "experimental", "mean_plddt": None}})  # no pockets
    d = dossier_node(st)["target_dossier"]["ORF9"]
    assert d["target_class"] is None                               # unknown, not invented
    assert d["pocket_descriptors"]["available"] is False
    assert d["pocket_descriptors"]["metal_ions_present"] is False
