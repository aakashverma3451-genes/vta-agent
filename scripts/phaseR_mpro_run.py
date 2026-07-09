"""Run SARS-CoV-2 Mpro (non-covalent) through the real Phase-R reasoning workflow.

Same pattern as scripts/phaseR_tilv_run.py, for the powered target: this is the one target
in the frozen benchmark with real measured inactives and informative (non-[0,1]) CIs.
"""
from __future__ import annotations

import json
from pathlib import Path

from vta.nodes.dossier import dossier_node
from vta.nodes.triage import triage_router_node
from vta.nodes.verification import build_verdict
from vta.state import new_state

OUT = Path("outputs/phaseR/mpro_run.json")


def run() -> dict:
    st = new_state("mpro-run", "mpro-run")
    st["structures"] = {"MPRO": {"method": "experimental", "source": "7L11:A", "mean_plddt": None}}
    st["pockets"] = {"MPRO": [{"id": 1, "center": [0, 0, 0], "volume": 500.0,
                               "druggability": 0.7, "conservation": 0.5}]}

    st = dossier_node(st)
    st = triage_router_node(st)
    triage = st["triage_decision"]["MPRO"]

    st["lead_candidates"] = [{"protein": "MPRO", "ligand": "probe", "dG": -8.0, "le": -0.3,
                              "conservation": 0.5, "score": 0.9}]
    verdict = build_verdict(st)

    result = {"target": "SARS-CoV-2 Mpro (7L11 chain A, non-covalent)",
              "triage": triage, "verification": verdict}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    r = run()
    print(f"Router decision: {r['triage']['decision']}")
    print(f"Verification verdict: {r['verification']['verdict']}")
    print(f"-> wrote {OUT}")
