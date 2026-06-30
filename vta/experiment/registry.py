"""Prospective prediction registry for locked experimental validation."""
from __future__ import annotations

from datetime import datetime, timezone

from vta.provenance import input_hash


def lock_predictions(leads: list[dict], run_id: str, top_n: int = 10) -> dict:
    locked = [
        {
            "rank": i + 1,
            "ligand": lead.get("ligand"),
            "ligand_id": lead.get("ligand_id"),
            "score": lead.get("score"),
            "input_hash": input_hash(lead),
        }
        for i, lead in enumerate(leads[:top_n])
    ]
    return {
        "run_id": run_id,
        "locked_at": datetime.now(timezone.utc).isoformat(),
        "status": "locked",
        "predictions": locked,
        "results": {},
    }


def ingest_result(registry: dict, ligand: str, result: dict) -> dict:
    if registry.get("status") != "locked":
        raise ValueError("registry must be locked before ingesting results")
    known = {p["ligand"] for p in registry.get("predictions") or []}
    if ligand not in known:
        raise ValueError(f"ligand was not prospectively locked: {ligand}")
    results = dict(registry.get("results") or {})
    if ligand in results:
        raise ValueError(f"experimental result already recorded for {ligand}")
    enriched = dict(result)
    if enriched.get("EC50") and enriched.get("CC50"):
        enriched["selectivity_index"] = round(float(enriched["CC50"]) / float(enriched["EC50"]), 3)
    enriched["ingested_at"] = datetime.now(timezone.utc).isoformat()
    results[ligand] = enriched
    return {**registry, "results": results}
