"""pockets_node (MOCK) — binding-pocket detection, FPocket-shaped output.

Plan Task 2.2. Mocked in Phase 1, real FPocket in Phase 2. The trick: the mock
returns records in EXACTLY the shape real FPocket produces — each pocket has a
druggability score (0–1), volume, geometric center, and a conservation score — so
the docking and ranking code downstream is built against the real shape and the
Phase-2 swap is invisible.

Determinism: pockets are generated from a seed derived from (run_id, protein), so a
given run is reproducible (the auditability thesis) while different proteins/runs
still differ. Real FPocket is of course deterministic per structure.
"""
from __future__ import annotations

import hashlib
import random

from vta.state import VTAState

_N_POCKETS = 3


def _seed(*parts: str) -> int:
    """Stable integer seed from string parts (hash() is salted per-process)."""
    h = hashlib.md5("|".join(parts).encode()).hexdigest()
    return int(h[:8], 16)


def _pockets_for(run_id: str, protein: str) -> list[dict]:
    rng = random.Random(_seed(run_id, protein))
    pockets = []
    for i in range(_N_POCKETS):
        pockets.append({
            "id": i + 1,
            "center": [round(rng.uniform(-15, 15), 1) for _ in range(3)],
            "volume": round(rng.uniform(350, 900), 0),
            "druggability": round(rng.uniform(0.30, 0.85), 2),
            "conservation": round(rng.uniform(0.30, 0.95), 2),
        })
    # FPocket ranks pockets best-first; sort by druggability descending.
    pockets.sort(key=lambda p: p["druggability"], reverse=True)
    for i, p in enumerate(pockets):
        p["id"] = i + 1
    return pockets


def pockets_node_mock(state: VTAState) -> VTAState:
    pockets = {}
    for name, struct in (state.get("structures") or {}).items():
        if not struct.get("pdb_path"):
            # no structure (refused / failed fold) → nothing to find pockets in
            state["audit_trail"].append(
                f"[MOCK] Pockets[{name}]: skipped (no structure)")
            continue
        plist = _pockets_for(state["run_id"], name)
        pockets[name] = plist
        state["audit_trail"].append(
            f"[MOCK] Pockets[{name}]: {len(plist)} found, "
            f"top druggability {plist[0]['druggability']}"
        )
    state["pockets"] = pockets
    state["versions"]["pockets"] = "MOCK-fpocket-shaped"
    return state
