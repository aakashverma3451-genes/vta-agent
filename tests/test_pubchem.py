"""PubChem resolver tests — hermetic (no network; the `_get` seam is monkeypatched)."""
from __future__ import annotations

import requests

import vta.data.pubchem as P


def _prop_table(cid, smiles):
    return {"PropertyTable": {"Properties": [{"CID": cid, "CanonicalSMILES": smiles}]}}


def test_resolve_smiles_returns_cid_and_smiles(monkeypatch):
    monkeypatch.setattr(P, "_get", lambda name, timeout=20: _prop_table(2244, "CC(=O)O"))
    hit = P.resolve_smiles("aspirin")
    assert hit == {"cid": 2244, "smiles": "CC(=O)O", "source": "PubChem"}


def test_resolve_smiles_404_is_definitive_miss(monkeypatch):
    # _get already maps 404 -> {}; resolve must return None without raising.
    monkeypatch.setattr(P, "_get", lambda name, timeout=20: {})
    assert P.resolve_smiles("not-a-real-compound") is None


def test_resolve_smiles_empty_name_returns_none():
    assert P.resolve_smiles("   ") is None


def test_resolve_smiles_transient_failure_returns_none(monkeypatch):
    monkeypatch.setattr(P, "_get",
                        lambda name, timeout=20: (_ for _ in ()).throw(
                            requests.RequestException("503")))
    monkeypatch.setattr(P.time, "sleep", lambda s: None)   # don't actually back off
    assert P.resolve_smiles("ribavirin") is None


def test_resolve_smiles_missing_smiles_field_returns_none(monkeypatch):
    monkeypatch.setattr(P, "_get",
                        lambda name, timeout=20: {"PropertyTable": {"Properties": [{"CID": 1}]}})
    assert P.resolve_smiles("x") is None


def test_ligands_uses_pubchem_when_chembl_misses(monkeypatch, tmp_path):
    """A ChEMBL miss is recovered from PubChem instead of being dropped."""
    import vta.data.ligands as L
    monkeypatch.setattr(L, "_fetch_one", lambda name, **kw: None)         # ChEMBL misses
    monkeypatch.setattr(L, "_resolve_via_pubchem",
                        lambda name: {"chembl_id": None, "pubchem_cid": 99,
                                      "name": name.title(), "smiles": "CCO",
                                      "max_phase": None, "positive_control": False,
                                      "source": "PubChem"})
    cache = tmp_path / "lig.json"
    out = L.refresh_ligands_from_chembl(names=["ETHANOL"], cache_path=str(cache), delay=0)
    assert len(out) == 1 and out[0]["source"] == "PubChem"
    assert out[0]["smiles"] == "CCO"
