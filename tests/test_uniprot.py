"""UniProt resolver tests — hermetic (no network; seams monkeypatched)."""
from __future__ import annotations

import vta.data.uniprot as U


def test_resolve_accession_returns_best_hit(monkeypatch, tmp_path):
    cache = tmp_path / "c.json"
    monkeypatch.setattr(U, "_search",
                        lambda q, timeout=20: [{"primaryAccession": "P03431"}])
    assert U.resolve_accession("PB1 influenza A", cache_path=str(cache)) == "P03431"


def test_resolve_accession_no_match_returns_none(monkeypatch, tmp_path):
    cache = tmp_path / "c.json"
    monkeypatch.setattr(U, "_search", lambda q, timeout=20: [])
    assert U.resolve_accession("nonsense xyz", cache_path=str(cache)) is None


def test_resolve_accession_is_cached_and_skips_network(monkeypatch, tmp_path):
    cache = tmp_path / "c.json"
    calls = {"n": 0}
    def counting_search(q, timeout=20):
        calls["n"] += 1
        return [{"primaryAccession": "Q99999"}]
    monkeypatch.setattr(U, "_search", counting_search)

    first = U.resolve_accession("PB2", cache_path=str(cache))
    second = U.resolve_accession("PB2", cache_path=str(cache))
    assert first == second == "Q99999"
    assert calls["n"] == 1                     # second call served from cache

    # A negative result is cached too — no repeat network for a hopeless query.
    monkeypatch.setattr(U, "_search", lambda q, timeout=20: [])
    assert U.resolve_accession("ghost", cache_path=str(cache)) is None
    monkeypatch.setattr(U, "_search",
                        lambda q, timeout=20: (_ for _ in ()).throw(AssertionError("hit net")))
    assert U.resolve_accession("ghost", cache_path=str(cache)) is None


def test_resolve_accession_empty_query_returns_none(tmp_path):
    assert U.resolve_accession("  ", cache_path=str(tmp_path / "c.json")) is None


def test_transient_failure_is_not_cached(monkeypatch, tmp_path):
    cache = tmp_path / "c.json"
    import requests
    monkeypatch.setattr(U, "_search",
                        lambda q, timeout=20: (_ for _ in ()).throw(
                            requests.RequestException("boom")))
    monkeypatch.setattr(U.time, "sleep", lambda s: None)   # don't actually back off
    assert U.resolve_accession("PA", cache_path=str(cache)) is None
    # transient miss must NOT be written to cache (so a later run can retry)
    assert not cache.exists() or "PA" not in U._load_cache(str(cache))


def test_fetch_sequence_strips_fasta_header(monkeypatch):
    monkeypatch.setattr(U, "_fasta", lambda acc, timeout=20: ">sp|P1|X\nMKT\nVLL\n")
    assert U.fetch_sequence("P1") == "MKTVLL"


def test_fetch_sequence_network_failure_returns_none(monkeypatch):
    import requests
    monkeypatch.setattr(U, "_fasta",
                        lambda acc, timeout=20: (_ for _ in ()).throw(
                            requests.RequestException("down")))
    assert U.fetch_sequence("P1") is None
