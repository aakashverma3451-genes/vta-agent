"""Tests for the HTML report generator (pure render + file write, proceed + defer)."""
from __future__ import annotations

from vta.report import render_report, write_report
from vta.state import new_state


def _proceed_state() -> dict:
    st = new_state("rep-proceed", "rep-proceed")
    st["taxon_result"] = {
        "genus": "Tilapinevirus", "species": "Tilapia lake virus", "strain_id": "TiLV",
        "confidence_basis": "aa identity, not calibrated",
        "taxonagent_version": "0.1.0", "kg_version": "ICTV-VMR-MSL40",
    }
    st["classification_confidence"] = 98.0
    st["route"] = "proceed"
    st["structures"] = {"PB1": {"method": "experimental", "source": "8PSO:B", "mean_plddt": None}}
    st["lead_candidates"] = [
        {"ligand": "Ribavirin", "ligand_id": "CHEMBL1643", "protein": "PB1",
         "dG": -6.79, "le": -0.399, "conservation": 0.5, "score": 0.6, "positive_control": True},
        {"ligand": "Lopinavir", "ligand_id": "CHEMBL729", "protein": "PB1",
         "dG": -8.8, "le": -0.191, "conservation": 0.5, "score": 0.498, "positive_control": False},
    ]
    st["audit_trail"] = ["TaxonAgent: ...", "Router: PROCEED", "Vina: 8 real docks"]
    st["versions"] = {"docking": "AutoDock Vina 1.2.5"}
    return st


def test_render_proceed_has_all_sections():
    html = render_report(_proceed_state())
    assert html.startswith("<!doctype html>")
    assert "Tilapia lake virus" in html
    assert "HIGH CONFIDENCE" in html and "98.0%" in html
    assert "experimental PDB 8PSO:B" in html
    assert "Ribavirin" in html and "★ control" in html      # control highlighted
    assert 'class="ctrl"' in html                            # control row styled
    assert "AutoDock Vina 1.2.5" in html
    assert "<details>" in html and "Router: PROCEED" in html  # collapsible audit


def test_render_flag_badge_is_amber():
    st = _proceed_state(); st["classification_confidence"] = 88.0; st["route"] = "flag"
    assert "FLAG CONFIDENCE" in render_report(st) and "badge a" in render_report(st)


def test_render_defer_has_no_leads():
    st = new_state("rep-defer", "rep-defer")
    st["taxon_result"] = {"genus": "X", "species": "Y", "confidence_basis": ""}
    st["classification_confidence"] = 40.0
    st["route"] = "defer"
    st["audit_trail"] = ["Router: DEFER"]
    html = render_report(st)
    assert "DEFERRED" in html and "badge r" in html
    assert "deferred to a human expert" in html               # empty-leads message


def test_write_report_creates_file(tmp_path):
    path = write_report(_proceed_state(), out_dir=str(tmp_path))
    assert path.endswith("rep-proceed_report.html")
    assert "<!doctype html>" in open(path).read()
