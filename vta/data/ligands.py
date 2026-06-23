"""vta.data.ligands — the real ligand library (Phase 2 drug-design data source).

Replaces the Phase-1 `MOCK_LIGANDS` placeholders with real compounds pulled from
**ChEMBL** (https://www.ebi.ac.uk/chembl/). Each ligand carries its ChEMBL id,
preferred name, canonical SMILES (needed for real docking ligand prep), and
`max_phase` (4 = approved drug). The fetched set is cached to a committed JSON so
the pipeline — and the tests — run offline; that cache file is VTA-Agent's first
real drug-design data asset.

Design notes:
  • We resolve a CURATED list of approved antivirals by name (the ChEMBL
    `pref_name__iexact` lookup, verified to return id + SMILES), rather than guessing
    ATC-filter syntax — every entry is reproducible and checkable.
  • The 5 nucleoside-analog RdRp inhibitors are flagged as POSITIVE CONTROLS: when
    real docking lands they must rank in the top decile against a polymerase target,
    or the scoring is mis-calibrated.
  • Network access is isolated in `refresh_ligands_from_chembl()`. Normal runs call
    `load_ligands()`, which reads the cache (no network).
"""
from __future__ import annotations

import json
import os
import time
from typing import Optional

import requests

_CHEMBL_MOLECULE = "https://www.ebi.ac.uk/chembl/api/data/molecule"
_CACHE_PATH = os.path.join(os.path.dirname(__file__), "ligands_chembl.json")

# Nucleoside/nucleotide-analog RdRp inhibitors — Phase-2 positive controls.
POSITIVE_CONTROLS = [
    "REMDESIVIR", "FAVIPIRAVIR", "RIBAVIRIN", "SOFOSBUVIR", "MOLNUPIRAVIR",
]

# Curated approved/known antivirals to dock against a viral polymerase. Polymerase-
# relevant nucleos(t)ide analogs are prioritised; a few other-mechanism antivirals
# are included as expected NEGATIVES (should rank low against an RdRp pocket).
CURATED_ANTIVIRALS = POSITIVE_CONTROLS + [
    # other polymerase-acting nucleos(t)ide analogs
    "TENOFOVIR", "LAMIVUDINE", "ENTECAVIR", "ACYCLOVIR", "GANCICLOVIR",
    "VALACYCLOVIR", "ZIDOVUDINE", "EMTRICITABINE", "ABACAVIR",
    # HCV NS5B / NS5A
    "DACLATASVIR", "LEDIPASVIR", "VELPATASVIR",
    # influenza (non-polymerase / cap-dependent) — mechanism contrast
    "OSELTAMIVIR", "ZANAMIVIR", "BALOXAVIR MARBOXIL", "PERAMIVIR",
    # protease / other classes — expected negatives for a polymerase pocket
    "NIRMATRELVIR", "LOPINAVIR", "RITONAVIR", "DOLUTEGRAVIR", "EFAVIRENZ",
]


def _fetch_one(name: str, timeout: int = 20, retries: int = 5) -> Optional[dict]:
    """Resolve one compound name → ligand dict via ChEMBL, or None if not found.

    ChEMBL throttles rapid sequential queries (429 / non-JSON error pages / read
    timeouts), so we retry transient failures with exponential backoff. A genuine
    "no such molecule" returns None immediately (no retry).
    """
    for attempt in range(retries):
        try:
            r = requests.get(
                _CHEMBL_MOLECULE,
                params={"pref_name__iexact": name, "format": "json"},
                timeout=timeout,
            )
            if r.status_code == 429 or r.status_code >= 500:
                raise requests.RequestException(f"transient HTTP {r.status_code}")
            r.raise_for_status()
            mols = r.json().get("molecules", [])           # JSONDecodeError → retry
            if not mols:
                return None                                # definitively not found
            m = mols[0]
            smiles = (m.get("molecule_structures") or {}).get("canonical_smiles")
            if not smiles:
                return None
            return {
                "chembl_id": m.get("molecule_chembl_id"),
                "name": (m.get("pref_name") or name).title(),
                "smiles": smiles,
                "max_phase": m.get("max_phase"),
                "positive_control": name in POSITIVE_CONTROLS,
            }
        except (requests.RequestException, ValueError):    # ValueError ⊇ JSONDecodeError
            if attempt == retries - 1:
                return None
            time.sleep(1.5 ** attempt)                     # 1, 1.5, 2.25, 3.4 s
    return None


def _resolve_via_pubchem(name: str) -> Optional[dict]:
    """PubChem fallback for a ChEMBL miss: recover SMILES so the compound isn't dropped.

    Thin seam over `vta.data.pubchem` (monkeypatched in tests). Returns a ligand dict
    shaped like a ChEMBL entry but sourced from PubChem, or None if PubChem misses too.
    """
    from vta.data.pubchem import resolve_smiles
    hit = resolve_smiles(name)
    if not hit:
        return None
    return {
        "chembl_id": None,
        "pubchem_cid": hit["cid"],
        "name": name.title(),
        "smiles": hit["smiles"],
        "max_phase": None,                               # unknown from PubChem alone
        "positive_control": name in POSITIVE_CONTROLS,
        "source": "PubChem",
    }


def refresh_ligands_from_chembl(
    names=None, cache_path: str = _CACHE_PATH, delay: float = 0.5
) -> list[dict]:
    """Fetch the curated set from ChEMBL (PubChem fallback) and write the cache. NETWORK.

    `delay` is a politeness pause between molecules to stay under ChEMBL's rate limit.
    A name ChEMBL can't resolve is retried against PubChem before being recorded as
    genuinely missing — so coverage gaps in one source don't silently drop a compound.
    """
    names = names or CURATED_ANTIVIRALS
    ligands, missing = [], []
    for name in names:
        lig = _fetch_one(name) or _resolve_via_pubchem(name)
        (ligands.append(lig) if lig else missing.append(name))
        time.sleep(delay)
    sources = sorted({lig.get("source", "ChEMBL") for lig in ligands})
    with open(cache_path, "w") as fh:
        json.dump({"source": "+".join(sources) or "ChEMBL",
                   "ligands": ligands, "missing": missing}, fh, indent=2)
    return ligands


# Small built-in fallback so an un-cached, offline import still yields *something*
# (clearly fake) rather than crashing. Real runs use the committed cache.
_MOCK_FALLBACK = [
    {"chembl_id": f"MOCK_{i:03d}", "name": f"mock_ligand_{i}", "smiles": None,
     "max_phase": None, "positive_control": False}
    for i in range(20)
]


def load_ligands(cache_path: str = _CACHE_PATH) -> list[dict]:
    """Return the real ligand library from cache; fall back to mock if absent."""
    if os.path.exists(cache_path):
        with open(cache_path) as fh:
            data = json.load(fh)
        if data.get("ligands"):
            return data["ligands"]
    return _MOCK_FALLBACK
