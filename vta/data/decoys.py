"""vta.data.decoys — decoy library for the enrichment benchmark (SPEC #5).

A retrospective enrichment benchmark needs DECOYS: compounds presumed inactive against
the target, mixed with the known actives so we can measure how well the score separates
them. The rigorous source is property-matched, topologically-dissimilar decoys per DUD-E
(Mysinger et al. 2012, doi:10.1021/jm300687e) — matched on MW/logP/HBD/HBA/rot-bonds/
charge but structurally different, so separation isn't trivial.

`fetch_decoys` is the injectable seam (monkeypatched in tests). Its real path would call
a DUD-E / LUDe / property-matched-ZINC generator; offline it returns a COMMITTED demo
decoy pool so the benchmark reproduces without network. No pool → None, and the driver
labels the skip honestly.

⚠ HONESTY: the committed pool is a DEMONSTRATION set (mechanism-distinct approved
antivirals — protease / NS5A / cap-endonuclease inhibitors that do NOT bind the RdRp NTP
site), NOT a DUD-E property-matched set. Real benchmarking must populate the DUD-E path;
DUD-E itself has documented analog bias (LIT-PCBA is the less-biased alternative, and a
2025 audit found leakage even there). The benchmark output states this caveat.
"""
from __future__ import annotations

import os
from typing import Optional

_DECOY_DIR = os.path.join(os.path.dirname(__file__), "decoys_cache")


def _read_smi(path: str) -> list[tuple[str, str]]:
    """Read a .smi (`SMILES  name` per line; '#'/blank ignored) → [(smiles, name)]."""
    out: list[tuple[str, str]] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            smiles = parts[0]
            name = parts[1].strip() if len(parts) > 1 else smiles
            out.append((smiles, name))
    return out


def load_decoy_pool() -> list[dict]:
    """All committed demo decoys: [{name, smiles, positive_control=False}]."""
    if not os.path.isdir(_DECOY_DIR):
        return []
    pool: list[dict] = []
    for fname in sorted(os.listdir(_DECOY_DIR)):
        if not fname.endswith(".smi"):
            continue
        for smiles, name in _read_smi(os.path.join(_DECOY_DIR, fname)):
            pool.append({"name": name, "smiles": smiles, "positive_control": False})
    return pool


def fetch_decoys(active_smiles, n: int = 50) -> Optional[list[str]]:
    """Decoy SMILES for the given active(s); None when no decoys are available.

    Real path (not wired here): DUD-E / LUDe / property-matched ZINC generation, ~`n`
    decoys per active. Offline: returns the committed demo pool (capped at `n` per
    active), excluding any that coincide with an active. None → labelled skip downstream.
    """
    actives = [active_smiles] if isinstance(active_smiles, str) else list(active_smiles)
    pool = [d["smiles"] for d in load_decoy_pool()]
    if not pool:
        return None
    active_set = set(actives)
    decoys = [s for s in pool if s not in active_set]
    cap = max(1, n) * max(1, len(actives))
    return decoys[:cap] if decoys else None
