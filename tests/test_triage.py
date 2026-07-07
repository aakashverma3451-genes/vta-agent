"""Tests for the TriageRouter node (Phase R / R2). Hermetic; drives via the dossier."""
from __future__ import annotations

from vta.nodes.dossier import dossier_node
from vta.nodes.triage import triage_router_node
from vta.state import new_state


def _routed(structures, pockets=None, screened=None):
    st = new_state("triage-t", "triage-t")
    st["structures"] = structures
    if pockets is not None:
        st["pockets"] = pockets
    if screened is not None:
        st["screened_ligand_class"] = screened
    st = dossier_node(st)
    st = triage_router_node(st)
    return st


def _dec(st, protein):
    return st["triage_decision"][protein]["decision"]


def test_mpro_routes_to_full_dock():
    st = _routed({"MPRO": {"method": "experimental", "source": "7L11:A", "mean_plddt": None}})
    assert _dec(st, "MPRO") == "full_dock"


def test_hcv_ni_routes_to_annotate_only_via_benchmarkability():
    # HCV NS5B is un_benchmarkable in the frozen benchmark → annotate_only (not a docking rank)
    st = _routed({"NS5B": {"method": "experimental", "source": "2XI3:A", "mean_plddt": None}})
    assert _dec(st, "NS5B") == "annotate_only"
    assert any("un_benchmarkable" in r for r in st["triage_decision"]["NS5B"]["rationale"])


def test_tilv_pb1_still_full_dock():
    # metal-dependent polymerase biology, but only UNDERPOWERED + no metal in pocket → dockable
    st = _routed({"PB1": {"method": "experimental", "source": "8PSO:B", "mean_plddt": None}})
    assert _dec(st, "PB1") == "full_dock"


def test_low_pocket_plddt_predicted_defers_or_refuses():
    defer = _routed({"GP120": {"method": "esmfold", "mean_plddt": 60.0,
                               "structure_qc": {"status": "flagged", "issues": ["low_plddt"]}}})
    assert _dec(defer, "GP120") == "defer"
    refuse = _routed({"GPX": {"method": "esmfold", "mean_plddt": 40.0,
                              "structure_qc": {"status": "flagged", "issues": ["low_plddt"]}}})
    assert _dec(refuse, "GPX") == "refuse"


def test_metal_in_pocket_downgrades_to_annotate():
    st = _routed({"MPRO": {"method": "experimental", "source": "7L11:A", "mean_plddt": None}},
                 {"MPRO": [{"id": 1, "center": [0, 0, 0], "metal_ions": ["MG"]}]})
    assert _dec(st, "MPRO") == "annotate_only"


def test_screened_nucleotide_ligand_downgrades_full_dock_target():
    st = _routed({"MPRO": {"method": "experimental", "source": "7L11:A", "mean_plddt": None}},
                 screened={"MPRO": "nucleotide"})
    assert _dec(st, "MPRO") == "annotate_only"


def test_playbook_stub_marks_no_precedent():
    st = _routed({"MPRO": {"method": "experimental", "source": "7L11:A", "mean_plddt": None}})
    assert st["playbook_prior"]["MPRO"]["status"] == "no_precedent"
