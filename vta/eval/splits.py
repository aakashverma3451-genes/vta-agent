"""Deterministic benchmark split helpers."""
from __future__ import annotations

import hashlib


def scaffold_key(smiles: str | None) -> str:
    """Return a stable scaffold-like key without requiring RDKit."""
    s = smiles or ""
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold
        mol = Chem.MolFromSmiles(s)
        if mol:
            scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
            return scaffold or s
    except Exception:
        pass
    return hashlib.sha1(s.encode()).hexdigest()[:12]


def scaffold_split(entries: list[dict], test_fraction: float = 0.2) -> dict[str, list[dict]]:
    keys = sorted({scaffold_key(e.get("smiles")) for e in entries})
    n_test = max(1, round(len(keys) * test_fraction)) if keys else 0
    test_keys = set(keys[:n_test])
    return {
        "train": [e for e in entries if scaffold_key(e.get("smiles")) not in test_keys],
        "test": [e for e in entries if scaffold_key(e.get("smiles")) in test_keys],
    }


def time_split(entries: list[dict], cutoff_year: int) -> dict[str, list[dict]]:
    return {
        "train": [e for e in entries if int(e.get("year") or 0) <= cutoff_year],
        "test": [e for e in entries if int(e.get("year") or 9999) > cutoff_year],
    }
