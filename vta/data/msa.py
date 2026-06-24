"""vta.data.msa — committed homolog MSA cache for conservation scoring (SPEC #1).

`conservation.fetch_homolog_msa` calls `load_msa_for(sequence)` to get a gapped homolog
alignment whose row 0 is the target sequence. Building such an alignment needs a homolog
search + MAFFT (network + an aligner), which is NOT reproducible in CI — so we BUILD it
offline (see `scripts/build_conservation_msa.py`) and COMMIT the result here as aligned
FASTA under `msa/`. At run time we just load the cache and match it to the query by exact
ungapped-row-0 identity, so the validation gate is fully reproducible offline.

No cache match → None, and the conservation node keeps the labelled 0.5 placeholder.
"""
from __future__ import annotations

import os
from typing import Optional

_MSA_DIR = os.path.join(os.path.dirname(__file__), "msa")


def _read_afa(path: str) -> list[str]:
    """Read an aligned FASTA into row strings, preserving order (row 0 = first record)."""
    rows: list[str] = []
    cur: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    rows.append("".join(cur))
                    cur = []
            else:
                cur.append(line.strip())
    if cur:
        rows.append("".join(cur))
    return rows


def _ungapped(row: str) -> str:
    return row.replace("-", "").replace(".", "")


def load_msa_for(sequence: str, taxon: Optional[str] = None) -> Optional[list[str]]:
    """Return committed MSA rows whose row-0 ungapped sequence equals `sequence`, or None.

    `taxon` is accepted for API symmetry with the live path but unused for cache lookup
    (the sequence identity is the precise key).
    """
    if not sequence or not os.path.isdir(_MSA_DIR):
        return None
    for fname in sorted(os.listdir(_MSA_DIR)):
        if not fname.endswith((".afa", ".aln", ".fasta")):
            continue
        rows = _read_afa(os.path.join(_MSA_DIR, fname))
        if rows and _ungapped(rows[0]) == sequence:
            return rows
    return None
