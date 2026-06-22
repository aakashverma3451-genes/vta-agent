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
# conservation isn't computed yet (needs an MSA); neutral so ranking stays balanced.
_CONSERVATION_PLACEHOLDER = 0.5


# ── real FPocket ─────────────────────────────────────────────────────────────
def _fpocket_bin() -> str | None:
    return os.environ.get("FPOCKET_BIN") or shutil.which("fpocket")


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
    """Real FPocket if available, else mock. Skips proteins with no structure."""
    fpocket = _fpocket_bin()
    pockets, used_real = {}, False
    for name, struct in (state.get("structures") or {}).items():
        pdb = struct.get("pdb_path")
        if not pdb:
            state["audit_trail"].append(f"Pockets[{name}]: skipped (no structure)")
            continue
        if fpocket:
            try:
                plist = _pockets_fpocket(name, pdb, fpocket)
                used_real = True
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
    return state


# back-compat alias (older imports / tests)
pockets_node_mock = pockets_node
