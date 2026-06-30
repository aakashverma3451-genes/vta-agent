"""Scientific gate: do the known RdRp-inhibitor controls rank high on a real screen?

Runs the REAL chain on TiLV PB1 — experimental 8PSO structure → real FPocket pocket
→ all cached ligands through real AutoDock Vina → composite ranking — then checks
whether the nucleoside-analog positive controls (remdesivir, ribavirin, sofosbuvir,
molnupiravir) land near the top. If they don't, the docking box / receptor prep is
suspect and must be debugged before trusting the pipeline.

Requires: FPOCKET_BIN, VINA_BIN, network (RCSB), and the rdkit/meeko prep stack.
Run from the vta-agent dir:
    FPOCKET_BIN=... VINA_BIN=... python scripts/validate_controls.py
"""
from __future__ import annotations

import json

from vta.nodes.docking import docking_node
from vta.nodes.conservation import conservation_node
from vta.nodes.conservation_contacts import conservation_contacts_node
from vta.nodes.pockets import pockets_node
from vta.nodes.rank import rank_node
from vta.nodes.structure import structure_node
from vta.state import new_state


def main() -> None:
    st = new_state("validation", "validation")
    # PB1 only (single segment); sequence unused on the experimental path.
    st["extracted_proteins"] = {"PB1": {"sequence": "M" * 500, "length_aa": 500, "plddt": None}}

    st = structure_node(st)      # real: fetch 8PSO chain B
    st = pockets_node(st)        # real: FPocket
    st = conservation_node(st)   # real: per-pocket JSD over committed homolog MSA (SPEC #1)
    st = docking_node(st)        # real: Vina over all cached ligands
    st = conservation_contacts_node(st)  # real: ligand-contact-weighted JSD (SPEC #4)
    st = rank_node(st)

    rows = sorted(st["docking_results"], key=lambda r: r["score"], reverse=True)
    print(f"\nengine: {st['versions'].get('docking')} | pockets: {st['versions'].get('pockets')}")
    print(f"ligands: {st['versions'].get('ligands')}\n")
    print(f"{'rank':<5}{'ligand':<22}{'dG':>8}{'LE':>9}{'cons':>7}{'score':>9}  control")
    print("-" * 70)
    for i, r in enumerate(rows, 1):
        print(f"{i:<5}{r['ligand']:<22}{r['dG']:>8}{r['le']:>9}{r['conservation']:>7}"
              f"{r['score']:>9}  {'★ CONTROL' if r['positive_control'] else ''}")

    # Scientific gate: controls in the top half (top-5 of 9).
    top_n = max(5, len(rows) // 2)
    top = rows[:top_n]
    ctrls_total = sum(1 for r in rows if r["positive_control"])
    ctrls_top = sum(1 for r in top if r["positive_control"])
    best_nonctrl = next((r for r in rows if not r["positive_control"]), None)
    print("-" * 70)
    print(f"controls in top {top_n}: {ctrls_top}/{ctrls_total}")
    verdict = "PASS" if ctrls_top >= max(1, ctrls_total - 1) else "FAIL — investigate box/prep"
    print(f"VERDICT: {verdict}")

    json.dump(
        {"rows": rows, "controls_in_top": ctrls_top, "controls_total": ctrls_total,
         "top_n": top_n, "verdict": verdict},
        open("outputs/validation_controls.json", "w"), indent=2,
    )
    print("wrote outputs/validation_controls.json")


if __name__ == "__main__":
    main()
