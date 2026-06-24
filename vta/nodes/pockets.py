"""pockets_node — binding-pocket detection.

Phase 2b: REAL FPocket when the binary is available, else a real-shaped MOCK.

`pockets_node` auto-detects fpocket (`FPOCKET_BIN` env var or `fpocket` on PATH).
If found it runs the real tool and parses its output; if not, it degrades to the
mock (same don't-crash discipline as the structure node) so CI without fpocket
still passes. Either way the output shape is identical, so docking/ranking downstream
don't care which produced it.

Real FPocket gives each pocket a druggability score, volume, and geometric center
(centroid of its alpha spheres). It does NOT give conservation — that needs a
multiple-sequence alignment — so `conservation` is set to a neutral placeholder and
flagged as a remaining Phase-2 task.
"""
from __future__ import annotations

import hashlib
import os
import random
import re
import shutil
import subprocess
from collections import defaultdict

from vta.state import VTAState

_N_POCKETS = 3
_TOP_KEEP = 3
# Neutral seed for conservation; the downstream `conservation` node (SPEC #1) overwrites
# this with a real per-pocket JSD score when a homolog MSA is available, else it stays 0.5.
_CONSERVATION_PLACEHOLDER = 0.5

# Known catalytic/active sites, taken from the BOUND SUBSTRATE in the experimental
# structure — the box center where the relevant inhibitors actually act. This is
# "experimental-first for pockets": a substrate-marked active site beats a blind
# FPocket scan, exactly as a solved structure beats a prediction. Coordinates are in
# the experimental structure's frame (same frame as the chain we extract).
#   PB1: centroid of the bound CTP nucleotide (8PSO chain F) — the +1 NTP site, where
#        nucleotide-analog RdRp inhibitors (remdesivir, sofosbuvir, …) are incorporated.
EXPERIMENTAL_ACTIVE_SITE = {
    "PB1": {"center": [135.56, 115.60, 123.66], "source": "8PSO:F (bound CTP, NTP site)"},
}


# ── real FPocket ─────────────────────────────────────────────────────────────
def _fpocket_bin() -> str | None:
    from vta.toolconfig import find_tool
    return find_tool("fpocket", "FPOCKET_BIN")


def _parse_info(info_path: str) -> dict[int, dict]:
    """Parse <base>_info.txt → {pocket_id: {druggability, volume, score}}."""
    pockets: dict[int, dict] = {}
    cur = None
    for line in open(info_path):
        m = re.match(r"\s*Pocket\s+(\d+)", line)
        if m:
            cur = int(m.group(1))
            pockets[cur] = {}
        elif cur is not None and ":" in line:
            key, _, val = line.partition(":")
            key, val = key.strip(), val.strip()
            try:
                num = float(val)
            except ValueError:
                continue
            if key == "Druggability Score":
                pockets[cur]["druggability"] = round(num, 3)
            elif key == "Volume":
                pockets[cur]["volume"] = round(num, 1)
            elif key == "Score":
                pockets[cur]["score"] = round(num, 3)
    return pockets


def _pocket_centers(out_pdb_path: str) -> dict[int, list]:
    """Centroid of each pocket's alpha spheres (STP HETATM records, grouped by id)."""
    acc = defaultdict(lambda: [0.0, 0.0, 0.0, 0])
    for line in open(out_pdb_path):
        if line.startswith("HETATM") and line[17:20].strip() == "STP":
            pid = int(line[22:26])
            a = acc[pid]
            a[0] += float(line[30:38]); a[1] += float(line[38:46]); a[2] += float(line[46:54])
            a[3] += 1
    return {pid: [round(a[0] / a[3], 3), round(a[1] / a[3], 3), round(a[2] / a[3], 3)]
            for pid, a in acc.items() if a[3]}


def _pockets_fpocket(name: str, pdb_path: str, fpocket: str) -> list[dict]:
    """Run real fpocket on a structure and return the top pockets (fpocket order)."""
    workdir = os.path.dirname(os.path.abspath(pdb_path))
    subprocess.run([fpocket, "-f", os.path.basename(pdb_path)],
                   cwd=workdir, capture_output=True, check=True, timeout=600)
    stem = os.path.splitext(os.path.basename(pdb_path))[0]
    out = os.path.join(workdir, f"{stem}_out")
    info = _parse_info(os.path.join(out, f"{stem}_info.txt"))
    centers = _pocket_centers(os.path.join(out, f"{stem}_out.pdb"))

    pockets = []
    for pid in sorted(info):                       # fpocket ranks best-first
        if pid not in centers:
            continue
        pockets.append({
            "id": pid,
            "center": centers[pid],
            "volume": info[pid].get("volume", 0.0),
            "druggability": info[pid].get("druggability", 0.0),
            "score": info[pid].get("score", 0.0),
            "conservation": _CONSERVATION_PLACEHOLDER,   # TODO: MSA-based, Phase 2
        })
    return pockets[:_TOP_KEEP]


# ── P2Rank consensus (§3.1) ──────────────────────────────────────────────────
def _p2rank_bin() -> str | None:
    from vta.toolconfig import find_tool
    return find_tool("prank", "P2RANK_BIN")


def _parse_p2rank(csv_path: str) -> list[dict]:
    """Parse P2Rank predictions CSV → [{rank, score, center}]."""
    results = []
    try:
        with open(csv_path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("name"):
                    continue
                parts = [p.strip() for p in line.split(",")]
                if len(parts) < 9:
                    continue
                results.append({
                    "rank": int(parts[1]),
                    "score": float(parts[2]),
                    "probability": float(parts[3]),
                    "center": [round(float(parts[6]), 3),
                               round(float(parts[7]), 3),
                               round(float(parts[8]), 3)],
                })
    except Exception:
        pass
    return sorted(results, key=lambda r: r["rank"])


def _pockets_p2rank(pdb_path: str, p2rank: str, workdir: str) -> list[dict]:
    """Run P2Rank on a PDB; return parsed pocket list (may be empty on failure)."""
    out = os.path.join(workdir, "p2rank_out")
    subprocess.run(
        [p2rank, "predict", "-f", pdb_path, "-o", out],
        capture_output=True, timeout=300,
    )
    stem = os.path.splitext(os.path.basename(pdb_path))[0]
    csv = os.path.join(out, f"{stem}.pdb_predictions.csv")
    if not os.path.exists(csv):
        return []
    return _parse_p2rank(csv)


def _dist(a: list, b: list) -> float:
    return sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5


def _add_p2rank_consensus(pockets: list[dict], p2rank_pockets: list[dict],
                          threshold: float = 4.0) -> None:
    """Annotate each FPocket pocket with consensus=True if P2Rank agrees (in-place)."""
    for p in pockets:
        c = p.get("center") or []
        if not c or not p2rank_pockets:
            p["consensus"] = False
            continue
        p["consensus"] = any(_dist(c, r["center"]) <= threshold
                             for r in p2rank_pockets[:3])


# ── mock fallback (real-shaped, deterministic) ───────────────────────────────
def _seed(*parts: str) -> int:
    return int(hashlib.md5("|".join(parts).encode()).hexdigest()[:8], 16)


def _pockets_mock(run_id: str, name: str) -> list[dict]:
    rng = random.Random(_seed(run_id, name))
    pockets = [{
        "id": i + 1,
        "center": [round(rng.uniform(-15, 15), 1) for _ in range(3)],
        "volume": round(rng.uniform(350, 900), 0),
        "druggability": round(rng.uniform(0.30, 0.85), 2),
        "conservation": round(rng.uniform(0.30, 0.95), 2),
    } for i in range(_N_POCKETS)]
    pockets.sort(key=lambda p: p["druggability"], reverse=True)
    for i, p in enumerate(pockets):
        p["id"] = i + 1
    return pockets


# ── node ─────────────────────────────────────────────────────────────────────
def pockets_node(state: VTAState) -> VTAState:
    """Real FPocket + P2Rank consensus if available, else mock.

    P2Rank runs as a parallel consensus check alongside FPocket: if both tools
    locate a pocket within 4 Å, that pocket gets consensus=True — a stronger
    prior for the docking box. Experimental active sites skip both detectors
    (they're ground truth) and get consensus=True by default.
    """
    fpocket = _fpocket_bin()
    p2rank = _p2rank_bin()
    pockets, used_real = {}, False
    for name, struct in (state.get("structures") or {}).items():
        pdb = struct.get("pdb_path")
        if not pdb:
            state["audit_trail"].append(f"Pockets[{name}]: skipped (no structure)")
            continue
        # Experimental-first: known active site, consensus implicit.
        if name in EXPERIMENTAL_ACTIVE_SITE:
            site = EXPERIMENTAL_ACTIVE_SITE[name]
            pockets[name] = [{
                "id": 1, "center": site["center"], "volume": None,
                "druggability": 1.0, "conservation": _CONSERVATION_PLACEHOLDER,
                "consensus": True, "detectors": ["experimental"],
                "method": "experimental_active_site", "source": site["source"],
            }]
            used_real = True
            state["audit_trail"].append(
                f"ActiveSite[{name}]: experimental catalytic site {site['source']} "
                f"@ {site['center']}")
            continue
        if fpocket:
            try:
                plist = _pockets_fpocket(name, pdb, fpocket)
                used_real = True
                # P2Rank consensus: annotate each FPocket pocket.
                detectors = ["fpocket"]
                if p2rank:
                    try:
                        pr = _pockets_p2rank(pdb, p2rank,
                                             os.path.dirname(os.path.abspath(pdb)))
                        _add_p2rank_consensus(plist, pr)
                        detectors.append("p2rank")
                        n_con = sum(1 for p in plist if p.get("consensus"))
                        state["audit_trail"].append(
                            f"P2Rank[{name}]: {len(pr)} pockets → "
                            f"{n_con}/{len(plist)} FPocket pockets have consensus")
                    except Exception as e:
                        state["audit_trail"].append(
                            f"P2Rank[{name}]: WARNING failed ({e}); skipping consensus")
                        for p in plist:
                            p["consensus"] = False
                else:
                    for p in plist:
                        p["consensus"] = False
                for p in plist:
                    p["detectors"] = detectors
                pockets[name] = plist
                state["audit_trail"].append(
                    f"FPocket[{name}]: {len(plist)} pockets, "
                    f"top druggability {plist[0]['druggability'] if plist else 'n/a'}")
                continue
            except Exception as e:
                state["audit_trail"].append(
                    f"FPocket[{name}]: WARNING failed ({e}); using mock")
        plist = _pockets_mock(state["run_id"], name)
        pockets[name] = plist
        state["audit_trail"].append(
            f"[MOCK] Pockets[{name}]: {len(plist)} found, "
            f"top druggability {plist[0]['druggability']}")
    state["pockets"] = pockets
    state["versions"]["pockets"] = "fpocket" if used_real else "MOCK-fpocket-shaped"
    if p2rank and used_real:
        state["versions"]["pockets"] += "+p2rank"
    return state


# back-compat alias (older imports / tests)
pockets_node_mock = pockets_node
