"""Phase S (power) — dock the undocked activity-cliff members to power the cliff benchmark.

The cliff benchmark is bounded by the DOCKED pool: of 1,515 Moonshot Mpro compounds with a
measured IC50, only 764 were docked (Job A), giving 17 scorable cliff pairs. Mining cliffs over
ALL measured compounds reveals a ceiling of ~1,193 pairs; docking the ~345 undocked compounds
that appear in those pairs unlocks them (a 70× increase → enough to test the +0.29 Vina−2D gap
at significance).

This docks exactly those undocked cliff-members into the SAME committed 7L11 receptor and box as
Job A (one consistent protocol), resumable + hang-capped (180 s/dock), writing ΔG to a Phase-S
cache. No new fetch, no re-fold — the receptor pdbqt is already on disk. A failed prep → None; a
hung dock → TIMEOUT; both labelled, never fabricated. Re-run `phaseS_activity_cliffs` afterwards.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List

from scripts.phase11_mpro_30to1 import _dock_one
from vta.eval.cliffs import find_cliff_pairs, pic50_from_um
from vta.nodes.docking import _prep_ligand

STRAT = Path("vta/data/mpro/mpro_stratified.json")
JOBA_CACHE = Path("outputs/phase11/.mpro_30to1_dgcache.json")
PHASES_CACHE = Path("outputs/phaseS/.cliff_dgcache.json")
RECEPTOR = "structures/mpro_30to1_MPRO_primary_rec.pdbqt"   # Job A receptor (already prepped)
CENTER = [-21.815, -4.216, -27.984]                          # 7L11 His41/Cys145 dyad box center
_LIG_CACHE = "structures/.ligand_pdbqt"
SIM_THRESHOLD, DPIC50_THRESHOLD = 0.7, 1.0


def _measured_compounds() -> List[dict]:
    s = json.loads(STRAT.read_text())
    cpds = []
    for key in ("non_covalent_actives", "inactives"):
        for r in s.get(key, []):
            m = r.get("measure") or {}
            p = pic50_from_um(m.get("ic50_um"))
            smi = r.get("canonical_smiles") or r.get("smiles")
            if p is not None and smi:
                cpds.append({"id": r["id"], "smiles": smi, "pic50": round(p, 4)})
    return cpds


def _load(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except Exception:
        return {}


def _undocked_cliff_members() -> List[dict]:
    """The exact set of compounds to dock: cliff-pair members with no committed ΔG yet."""
    cpds = _measured_compounds()
    docked = {k for k, v in {**_load(JOBA_CACHE), **_load(PHASES_CACHE)}.items()
              if isinstance(v, (int, float))}
    pairs, _ = find_cliff_pairs(cpds, sim_threshold=SIM_THRESHOLD, dpic50_threshold=DPIC50_THRESHOLD)
    need_ids = set()
    for (i, j, _s, _d) in pairs:
        for k in (i, j):
            if cpds[k]["id"] not in docked:
                need_ids.add(cpds[k]["id"])
    by_id = {c["id"]: c for c in cpds}
    return [by_id[i] for i in sorted(need_ids)]


def run() -> dict:
    todo = _undocked_cliff_members()
    cache = _load(PHASES_CACHE)
    PHASES_CACHE.parent.mkdir(parents=True, exist_ok=True)
    done = fail = timeout = 0
    for n, c in enumerate(todo, 1):
        cid = c["id"]
        if cid in cache:
            continue
        lig = {"smiles": c["smiles"], "dock_cache_id": cid, "chembl_id": cid, "name": cid}
        try:
            lp = _prep_ligand(lig, _LIG_CACHE)
            if not lp:
                cache[cid] = None
                fail += 1
            else:
                dG = _dock_one(RECEPTOR, lp, CENTER, f"structures/phaseS_{cid}_pose.pdbqt")
                cache[cid] = dG
                if dG is None:
                    fail += 1
                else:
                    done += 1
        except Exception:
            cache[cid] = "TIMEOUT"
            timeout += 1
        PHASES_CACHE.write_text(json.dumps(cache))
        if n % 10 == 0:
            print(f"progress {n}/{len(todo)} docked={done} fail={fail} timeout={timeout}", flush=True)
    numeric = sum(isinstance(v, (int, float)) for v in cache.values())
    print(f"DONE: {len(todo)} targets, {numeric} numeric ΔG in Phase-S cache "
          f"(fail={fail} timeout={timeout})", flush=True)
    return {"n_todo": len(todo), "numeric_dG": numeric, "fail": fail, "timeout": timeout}


if __name__ == "__main__":
    run()
