"""Hermetic tests for Phase 9B-3 HCV NS5B NI decoy matching and provenance.

The ChEMBL network seam (``fetch_fn``) is replaced by an in-memory fake, so these
tests never touch the network.
"""
from __future__ import annotations

import pytest

from scripts.fetch_phase9_decoys import (
    PUBLICATION_GRADE_PER_ACTIVE,
    build_decoys,
    compute_properties,
    largest_fragment_smiles,
    property_match,
    summarize_matching_quality,
)

# A real NS5B NI active-site active (sofosbuvir triphosphate active form, salt-free).
ACTIVE_SMILES = "C[C@@]1(O)[C@H](O)[C@@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O[C@H]1n1ccc(=O)[nH]c1=O"
# Same active but written with triethylamine counterions, as ChEMBL sometimes exports.
ACTIVE_WITH_SALT = ACTIVE_SMILES + ".CCN(CC)CC.CCN(CC)CC"


def _active(idx: str, smiles: str = ACTIVE_SMILES) -> dict:
    return {"molecule_chembl_id": idx, "canonical_smiles": smiles, "pref_name": idx}


def test_largest_fragment_strips_salt():
    assert largest_fragment_smiles(ACTIVE_WITH_SALT) == largest_fragment_smiles(ACTIVE_SMILES)
    assert "CCN(CC)CC" not in (largest_fragment_smiles(ACTIVE_WITH_SALT) or "")


def test_properties_ignore_counterions():
    bare = compute_properties(ACTIVE_SMILES)
    salted = compute_properties(ACTIVE_WITH_SALT)
    assert bare is not None and salted is not None
    # Salt stripping means identical descriptors despite the counterions.
    assert bare["mw"] == salted["mw"]
    assert bare["charge"] == salted["charge"]


def test_property_match_window():
    props = compute_properties(ACTIVE_SMILES)
    near = dict(props)
    near.update({"mw": props["mw"] + 10, "alogp": props["alogp"] + 0.5})
    assert property_match(props, near)
    far = dict(props)
    far["mw"] = props["mw"] + 500
    assert not property_match(props, far)


def _fake_fetch_factory(molecules: list[dict]):
    """Return a fetch_fn that always yields the same candidate molecules."""
    def fetch_fn(params: dict) -> list[dict]:
        return [
            {
                "molecule_chembl_id": m["id"],
                "pref_name": m.get("name"),
                "molecule_structures": {"canonical_smiles": m["smiles"]},
            }
            for m in molecules
        ]
    return fetch_fn


def test_build_decoys_matches_and_excludes_active_scaffold():
    actives = [_active("ACT1")]
    active_scaffold_smiles = largest_fragment_smiles(ACTIVE_SMILES)
    candidates = [
        # Scaffold-identical to the active -> must be excluded even if property-matched.
        {"id": "DECOY_SAME_SCAFFOLD", "smiles": ACTIVE_SMILES, "name": "same"},
        # Property-matched, scaffold-distinct triphosphate -> should be kept.
        {"id": "DECOY_OK", "name": "ok",
         "smiles": "Nc1nc(=O)n(cc1)[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O"},
        # Way off on MW -> filtered out by property match.
        {"id": "DECOY_FAR", "smiles": "CCO", "name": "far"},
    ]
    prov = build_decoys(actives, per_active=30, pages=1,
                        fetch_fn=_fake_fetch_factory(candidates), delay=0.0)
    ids = {d["chembl_id"] for d in prov["decoys"]}
    assert "DECOY_OK" in ids
    assert "DECOY_SAME_SCAFFOLD" not in ids  # scaffold-distinct enforced
    assert "DECOY_FAR" not in ids            # property window enforced
    assert active_scaffold_smiles is not None


def test_build_decoys_dedups_and_skips_active_ids():
    actives = [_active("ACT1"), _active("ACT2")]
    candidates = [
        {"id": "ACT1", "smiles": ACTIVE_SMILES, "name": "is-active"},  # an active id
        {"id": "DECOY_OK", "name": "ok",
         "smiles": "Nc1nc(=O)n(cc1)[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O"},
    ]
    prov = build_decoys(actives, per_active=30, pages=1,
                        fetch_fn=_fake_fetch_factory(candidates), delay=0.0)
    ids = [d["chembl_id"] for d in prov["decoys"]]
    assert "ACT1" not in ids
    assert ids.count("DECOY_OK") == 1  # deduped across the two actives


def test_underpowered_set_labelled_demonstration_grade():
    actives = [_active("ACT1")]
    candidates = [
        {"id": "DECOY_OK", "name": "ok",
         "smiles": "Nc1nc(=O)n(cc1)[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O"},
    ]
    prov = build_decoys(actives, per_active=30, pages=1,
                        fetch_fn=_fake_fetch_factory(candidates), delay=0.0)
    assert prov["n_decoys"] < PUBLICATION_GRADE_PER_ACTIVE
    assert prov["publication_grade"] is False
    assert prov["grade"] == "demonstration-grade"


def test_network_failure_degrades_with_provenance():
    def failing_fetch(params: dict):
        import requests
        raise requests.RequestException("simulated outage")

    prov = build_decoys([_active("ACT1")], per_active=30, pages=1,
                        fetch_fn=failing_fetch, delay=0.0)
    assert prov["n_decoys"] == 0
    assert prov["n_failures"] == 1
    assert "simulated outage" in prov["failures"][0]["error"]
    assert prov["publication_grade"] is False


def test_matching_quality_report_shape():
    actives = [_active("ACT1")]
    candidates = [
        {"id": "DECOY_OK", "name": "ok",
         "smiles": "Nc1nc(=O)n(cc1)[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O"},
    ]
    prov = build_decoys(actives, per_active=30, pages=1,
                        fetch_fn=_fake_fetch_factory(candidates), delay=0.0)
    quality = summarize_matching_quality(prov, actives)
    assert quality["n_actives"] == 1
    assert quality["scaffold_distinct_from_actives"] is True
    for prop in ("mw", "alogp", "hba", "hbd", "rtb", "charge"):
        assert prop in quality["property_distribution"]
        assert "actives" in quality["property_distribution"][prop]
        assert "decoys" in quality["property_distribution"][prop]
    assert quality["grade"] == prov["grade"]
