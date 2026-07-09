"""External database adapter tests — hermetic, no network."""
from __future__ import annotations

from vta.data import databases as db
from vta.data import external_databases as ext


DEDICATED_KEYS = {
    "rcsb_pdb",
    "esm_atlas",
    "alphafold",
    "uniprot",
    "homolog_msa",
    "chembl",
    "pubchem",
}


def test_every_non_dedicated_database_has_adapter():
    external_keys = {d.key for d in db.DATABASES if d.key not in DEDICATED_KEYS}
    assert external_keys == set(ext.ADAPTERS)


def test_adapter_urls_are_deterministic():
    assert ext.build_url("pdbe", "8PSO").endswith("/8PSO")
    assert "query.term=ribavirin" in ext.build_url("clinicaltrials", "ribavirin")
    assert ext.build_url("swissadme") == "http://www.swissadme.ch"


def test_json_fetch_skips_non_json_adapters(monkeypatch):
    def fail_get(*_args, **_kwargs):
        raise AssertionError("non-json adapter should not call requests.get")

    monkeypatch.setattr(ext.requests, "get", fail_get)
    assert ext.fetch_json("swissadme", "CCO") is None


def test_json_fetch_uses_adapter_url(monkeypatch):
    calls = []

    class Response:
        status_code = 200

        def raise_for_status(self):
            return None

        def json(self):
            return {"ok": True}

    def fake_get(url, timeout):
        calls.append((url, timeout))
        return Response()

    monkeypatch.setattr(ext.requests, "get", fake_get)
    assert ext.fetch_json("pdbe", "8PSO", timeout=7) == {"ok": True}
    assert calls == [("https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/8PSO", 7)]


def test_adapter_summary_is_serializable_shape():
    summary = ext.adapter_summary()
    assert summary
    assert {"key", "name", "access", "base_url"} <= set(summary[0])
