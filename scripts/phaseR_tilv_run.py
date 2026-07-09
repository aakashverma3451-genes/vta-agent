"""Run TiLV PB1 through the real Phase-R reasoning workflow (dossier -> triage -> verification).

This is not a new docking run — it exercises the shipped R1/R2/R4 node functions against the
committed TiLV PB1 benchmark artifacts (outputs/phase8/*, outputs/phase11/baselines_and_gate.json,
outputs/phase11/redock_validation.json) to show what the new architecture concludes about this
target, as opposed to the pre-Phase-R pipeline which would dock and report a ranking with no
gate on the claim.
"""
from __future__ import annotations

import json
from pathlib import Path

from vta.nodes.dossier import dossier_node
from vta.nodes.triage import triage_router_node
from vta.nodes.verification import build_verdict
from vta.state import new_state

OUT = Path("outputs/phaseR/tilv_pb1_run.json")


def run() -> dict:
    st = new_state("tilv-run", "tilv-run")
    st["structures"] = {"PB1": {"method": "experimental", "source": "8PSO:B", "mean_plddt": None}}
    st["pockets"] = {"PB1": [{"id": 1, "center": [0, 0, 0], "volume": 500.0,
                              "druggability": 0.7, "conservation": 0.5}]}

    st = dossier_node(st)
    st = triage_router_node(st)
    triage = st["triage_decision"]["PB1"]

    st["lead_candidates"] = [{"protein": "PB1", "ligand": "Ribavirin", "dG": -7.5, "le": -0.3,
                              "conservation": 0.5, "score": 0.635}]
    verdict = build_verdict(st)

    result = {"target": "TiLV PB1 (8PSO chain B)", "triage": triage, "verification": verdict}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    r = run()
    print(f"Router decision: {r['triage']['decision']}")
    print(f"Verification verdict: {r['verification']['verdict']}")
    print(f"-> wrote {OUT}")
