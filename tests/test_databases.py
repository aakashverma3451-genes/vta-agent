"""Drug-design database registry tests — hermetic (no network).

Guards the catalog in `vta/data/databases.py` against drift: unique keys, honest
integration flags, and consistency with the modules that actually pull data.
"""
from __future__ import annotations

import importlib

import pytest

from vta.data import databases as db
from vta.data.databases import Category, Status


def test_keys_are_unique():
    keys = [d.key for d in db.DATABASES]
    assert len(keys) == len(set(keys)), "duplicate database key in catalog"


def test_every_category_is_represented():
    present = {d.category for d in db.DATABASES}
    assert present == set(Category), "some Category has no databases listed"


def test_lookup_roundtrips():
    for d in db.DATABASES:
        assert db.get(d.key) is d
    with pytest.raises(KeyError):
        db.get("no_such_database")


def test_integrated_sources_name_a_real_importable_module():
    integrated = db.integrated()
    assert integrated, "expected at least the structure + ligand sources wired in"
    for d in integrated:
        # INTEGRATED must declare WHERE it's wired, and that module must import.
        assert d.wired_in, f"{d.key} is INTEGRATED but names no module"
        importlib.import_module(d.wired_in)


def test_only_integrated_entries_claim_a_module():
    # Honesty guard: don't claim a wiring module unless status is INTEGRATED.
    for d in db.DATABASES:
        if d.wired_in:
            assert d.status is Status.INTEGRATED, (
                f"{d.key} names a module but isn't marked INTEGRATED")


def test_known_real_sources_are_integrated():
    # These two are genuinely pulled from today; lock that in.
    assert db.get("rcsb_pdb").status is Status.INTEGRATED
    assert db.get("chembl").status is Status.INTEGRATED
    assert db.get("chembl").wired_in == "vta.data.ligands"
    assert db.get("rcsb_pdb").wired_in == "vta.nodes.structure"


def test_summary_mentions_integrated_sources():
    text = db.summary()
    assert "ChEMBL" in text and "PDB (RCSB)" in text
    assert "✅" in text
