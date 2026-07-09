"""Small provenance helpers for numeric VTA outputs."""
from __future__ import annotations

import hashlib
import json
from typing import Any


def input_hash(value: Any) -> str:
    """Stable short hash for JSON-like input data."""
    payload = json.dumps(value, sort_keys=True, default=str, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def score_provenance(tool: str, version: str | None = None, *,
                     seed: int | str | None = None, inputs: Any = None) -> dict:
    """Provenance block attached next to score-like fields."""
    out = {"tool": tool}
    if version:
        out["version"] = version
    if seed is not None:
        out["seed"] = seed
    if inputs is not None:
        out["input_hash"] = input_hash(inputs)
    return out


def run_provenance(state: dict) -> dict:
    """Machine-readable run provenance suitable for RO-Crate/PROV-O export."""
    return {
        "run_id": state.get("run_id"),
        "versions": state.get("versions") or {},
        "input_hash": input_hash({
            "genome_fasta": state.get("genome_fasta"),
            "taxon_result": state.get("taxon_result"),
            "structures": state.get("structures"),
        }),
        "audit_hash": input_hash(state.get("audit_trail") or []),
        "provenance_model": "VTA lightweight PROV-O/RO-Crate compatible block",
    }
