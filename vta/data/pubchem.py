"""vta.data.pubchem — resolve compounds via PubChem PUG-REST (ligand-source fallback).

ChEMBL is the primary ligand library (`vta.data.ligands`), but it doesn't have every
compound under the exact name we query. PubChem is the broadest public small-molecule
database and exposes a clean REST API (PUG-REST), so it's the natural fallback: when a
curated antiviral name misses in ChEMBL, resolve its canonical SMILES here instead of
dropping it.

Mirrors the discipline of `vta.data.ligands` / `vta.data.uniprot`:
  • Network access is isolated in an injectable seam (`_get`) so tests stay offline.
  • A genuine "no such compound" returns None; transient failures (429/5xx/timeout)
    retry with backoff.

Catalogued in `vta.data.databases` as key "pubchem".
"""
from __future__ import annotations

import time
from typing import Optional

import requests

# PUG-REST: name -> property table (canonical SMILES + CID).
_PUG_REST = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/"
    "property/CanonicalSMILES/JSON"
)


def _get(name: str, timeout: int = 20) -> dict:
    """Fetch the PUG-REST property JSON for a compound name. Raises on HTTP error.

    A 404 (PUG's "no record" for an unknown name) is surfaced as a definitive miss by
    the caller, NOT retried.
    """
    r = requests.get(_PUG_REST.format(name=requests.utils.quote(name)), timeout=timeout)
    if r.status_code == 404:
        return {}                                        # definitive: no such compound
    if r.status_code == 429 or r.status_code >= 500:
        raise requests.RequestException(f"transient HTTP {r.status_code}")
    r.raise_for_status()
    return r.json()


def resolve_smiles(name: str, *, retries: int = 4) -> Optional[dict]:
    """Resolve a compound name to {cid, smiles, source} via PubChem, or None.

    None means "not found" (or unreachable after retries); the caller then drops the
    compound, exactly as it would for a ChEMBL miss.
    """
    name = (name or "").strip()
    if not name:
        return None

    for attempt in range(retries):
        try:
            props = (_get(name).get("PropertyTable", {}).get("Properties") or [])
            if not props:
                return None                              # definitive miss (incl. 404)
            p = props[0]
            smiles = p.get("CanonicalSMILES")
            if not smiles:
                return None
            return {"cid": p.get("CID"), "smiles": smiles, "source": "PubChem"}
        except (requests.RequestException, ValueError):  # ValueError ⊇ JSONDecodeError
            if attempt == retries - 1:
                return None
            time.sleep(1.5 ** attempt)                   # 1, 1.5, 2.25 s
    return None
