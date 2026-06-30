"""Hermetic tests for the Mpro covalent-warhead classifier (no network)."""
from __future__ import annotations

from vta.chem.warheads import classify_binding_mode, matched_warheads

# A real Moonshot nitrile fragment (row 0 of covid_submissions_all_info.csv), flagged covalent.
NITRILE = "N#Cc1ccccc1NC(=O)Cc1c[nH]c2ncccc12"
# A non-covalent amide-type Moonshot lead chemotype (no electrophilic warhead).
NONCOVALENT = "Cc1ccncc1NC(=O)Cc1cccc(Cl)c1"
ALDEHYDE = "O=CC(Cc1ccccc1)NC(=O)C1CC1"
KETOAMIDE = "O=C(N)C(=O)C(Cc1ccccc1)NC(=O)C1CC1"
CHLOROACETAMIDE = "O=C(CCl)Nc1ccccc1"


def test_substructure_detects_covalent_warheads():
    assert "nitrile" in matched_warheads(NITRILE)
    assert "aldehyde" in matched_warheads(ALDEHYDE)
    assert "alpha_ketoamide" in matched_warheads(KETOAMIDE)
    assert "haloacetamide" in matched_warheads(CHLOROACETAMIDE)


def test_noncovalent_lead_has_no_warhead():
    assert matched_warheads(NONCOVALENT) == []
    out = classify_binding_mode(NONCOVALENT)
    assert out["mode"] == "non_covalent"
    assert out["basis"] == "substructure"


def test_substructure_classifies_nitrile_as_covalent_without_annotation():
    out = classify_binding_mode(NITRILE)
    assert out["mode"] == "covalent"
    assert out["basis"] == "substructure"
    assert "nitrile" in out["warheads"]


def test_explicit_annotation_wins_over_substructure():
    # Moonshot explicitly flags this nitrile covalent.
    out = classify_binding_mode(NITRILE, annotation="True")
    assert out["mode"] == "covalent"
    assert out["basis"] == "annotation"
    # An authoritative non-covalent annotation is trusted even if a warhead is present
    # (e.g. an aryl nitrile that is a non-reactive substituent here).
    out2 = classify_binding_mode(NITRILE, annotation="False")
    assert out2["mode"] == "non_covalent"
    assert out2["basis"] == "annotation"
    assert "nitrile" in out2["warheads"]  # tension surfaced, not hidden


def test_blank_annotation_falls_through_to_substructure():
    out = classify_binding_mode(NONCOVALENT, annotation="")
    assert out["mode"] == "non_covalent"
    assert out["basis"] == "substructure"


def test_unparseable_is_ambiguous():
    out = classify_binding_mode("not a smiles %%%")
    assert out["mode"] == "ambiguous"
    assert out["basis"] == "unparseable"
