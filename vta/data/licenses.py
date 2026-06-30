"""Database license/access metadata registry."""
from __future__ import annotations

DATABASES = {
    "ChEMBL": {"license": "CC BY-SA 3.0", "access": "open"},
    "PDB": {"license": "CC0-like data policy", "access": "open"},
    "DrugBank": {"license": "restricted", "access": "explicit_opt_in"},
    "ZINC": {"license": "source-specific", "access": "open_with_terms"},
    "Enamine REAL": {"license": "commercial terms", "access": "explicit_opt_in"},
    "PrimeKG": {"license": "source-specific", "access": "open_with_terms"},
    "DRKG": {"license": "source-specific", "access": "open_with_terms"},
}


def validate_registry(registry: dict = DATABASES) -> tuple[bool, list[str]]:
    missing = [name for name, meta in registry.items() if not meta.get("license") or not meta.get("access")]
    return not missing, missing
