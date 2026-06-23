"""vta.data.uniprot — resolve UniProt accessions + sequences (drug-design data source).

This is the missing upstream half of the AlphaFold-DB structure fallback: AlphaFold DB
is keyed by UniProt accession, but nothing in the pipeline assigns one. Given a free-text
query (protein name + organism, e.g. "PB1 influenza A polymerase"), `resolve_accession`
asks the UniProt REST API for the best-matching accession, which the structure node then
feeds to AlphaFold DB. `fetch_fasta` returns the canonical sequence for an accession.

Mirrors the discipline of `vta.data.ligands`:
  • Network access is isolated in injectable seams (`_search`, `_fasta`) so tests stay
    offline by monkeypatching them.
  • Resolved accessions are memo-cached to a committed JSON so repeat runs — and CI —
    don't re-hit the API.
  • Transient failures (429 / 5xx / timeouts) retry with backoff; a genuine "no match"
    returns None immediately.

Catalogued in `vta.data.databases` as key "uniprot".
"""
from __future__ import annotations

import json
import os
import time
from typing import Optional

import requests

_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"
_CACHE_PATH = os.path.join(os.path.dirname(__file__), "uniprot_cache.json")


# ── injectable network seams (monkeypatched in tests) ────────────────────────
def _search(query: str, timeout: int = 20) -> list[dict]:
    """Query the UniProt search API; return the `results` list. Raises on HTTP error."""
    r = requests.get(
        _SEARCH_URL,
        params={
            "query": query,
            "format": "json",
            "size": 1,                                   # best hit only
            "fields": "accession,protein_name,organism_name,length",
        },
        timeout=timeout,
    )
    if r.status_code == 429 or r.status_code >= 500:
        raise requests.RequestException(f"transient HTTP {r.status_code}")
    r.raise_for_status()
    return r.json().get("results", [])


def _fasta(acc: str, timeout: int = 20) -> str:
    """Fetch the canonical FASTA for an accession. Raises on HTTP error."""
    r = requests.get(_ENTRY_URL.format(acc=acc), timeout=timeout)
    r.raise_for_status()
    return r.text


# ── cache helpers ────────────────────────────────────────────────────────────
def _load_cache(path: str) -> dict:
    if os.path.exists(path):
        try:
            with open(path) as fh:
                return json.load(fh)
        except (OSError, ValueError):
            return {}
    return {}


def _save_cache(cache: dict, path: str) -> None:
    with open(path, "w") as fh:
        json.dump(cache, fh, indent=2, sort_keys=True)


# ── public API ───────────────────────────────────────────────────────────────
def resolve_accession(
    query: str, *, cache_path: str = _CACHE_PATH, retries: int = 4
) -> Optional[str]:
    """Resolve a free-text protein query to its best UniProt accession, or None.

    Result is memo-cached by query string so repeat lookups don't hit the network.
    A definitive "no match" is cached as None too (so we don't retry a hopeless query
    every run). Transient HTTP failures retry with exponential backoff and are NOT
    cached.
    """
    query = (query or "").strip()
    if not query:
        return None

    cache = _load_cache(cache_path)
    if query in cache:                                   # includes cached None
        return cache[query]

    for attempt in range(retries):
        try:
            results = _search(query)
            acc = results[0].get("primaryAccession") if results else None
            cache[query] = acc                           # cache hit OR definitive miss
            _save_cache(cache, cache_path)
            return acc
        except (requests.RequestException, ValueError):
            if attempt == retries - 1:
                return None                              # transient → don't poison cache
            time.sleep(1.5 ** attempt)                   # 1, 1.5, 2.25 s
    return None


def fetch_sequence(acc: str) -> Optional[str]:
    """Return the canonical amino-acid sequence for an accession, or None on failure."""
    try:
        fasta = _fasta(acc)
    except (requests.RequestException, ValueError):
        return None
    seq = "".join(
        line.strip() for line in fasta.splitlines() if line and not line.startswith(">")
    )
    return seq or None
