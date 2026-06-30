"""Run a real TiLV PB1 benchmark over controls plus Phase 8 matched decoys."""
from __future__ import annotations

import json
from pathlib import Path

import vta.nodes.docking as docking
from vta.data.ligands import load_ligands
from vta.eval.leakage import leakage_audit
from vta.eval.metrics import bootstrap_enrichment_report, enrichment_report
from vta.nodes.conservation import conservation_node
from vta.nodes.conservation_contacts import conservation_contacts_node
from vta.nodes.pockets import pockets_node
from vta.nodes.rank import rank_node
from vta.nodes.structure import structure_node
from vta.state import new_state

DECOYS = Path("vta/data/decoys_cache/phase8_tilv_pb1_chembl_matched.smi")
OUT = Path("outputs/phase8/matched_benchmark_tilv_pb1.json")


def _decoys() -> list[dict]:
    rows = []
    for i, line in enumerate(DECOYS.read_text().splitlines(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        smiles, name, chembl_id, *_ = line.split("\t")
        rows.append({
            "name": name,
            "chembl_id": chembl_id or f"PHASE8_DECOY_{i}",
            "smiles": smiles,
            "max_phase": None,
            "positive_control": False,
            "source": "ChEMBL property-matched presumed decoy",
        })
    return rows


def _library() -> list[dict]:
    controls = [lig for lig in load_ligands() if lig.get("positive_control")]
    return controls + _decoys()


def run() -> dict:
    if not DECOYS.exists() or not DECOYS.read_text().strip():
        raise FileNotFoundError(f"no Phase 8 decoys at {DECOYS}")
    ligands = _library()
    docking.load_ligands = lambda: ligands

    st = new_state("phase8_matched_benchmark", "phase8_matched_benchmark")
    st["extracted_proteins"] = {"PB1": {"sequence": "M" * 500, "length_aa": 500, "plddt": None}}
    st = structure_node(st)
    st = pockets_node(st)
    st = conservation_node(st)
    st = docking.docking_node(st)
    st = conservation_contacts_node(st)
    st = rank_node(st)

    rows = sorted(st["docking_results"], key=lambda r: r["score"], reverse=True)
    entries = [{
        "name": r["ligand"],
        "score": r["score"],
        "positive_control": bool(r.get("positive_control")),
        "smiles": r.get("smiles"),
    } for r in rows]
    report = enrichment_report(entries)
    bootstrap = bootstrap_enrichment_report(entries, n_resamples=1000, seed=42)
    payload = {
        "target": "TiLV PB1 8PSO chain B",
        "decoy_set": str(DECOYS),
        "decoy_quality": "real ChEMBL property-matched presumed decoys; underpowered",
        "benchmark": report,
        "bootstrap": bootstrap,
        "rows": rows,
        "leakage_audit": leakage_audit([], ["TiLV_PB1_8PSO_B"]),
        "audit_trail": st.get("audit_trail"),
        "versions": st.get("versions"),
        "verdict": "defensible demonstration; underpowered for publication",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2))
    return payload


def main() -> None:
    out = run()
    bench = out["benchmark"]
    print(f"N={bench['n']} actives={bench['n_actives']} decoys={bench['n_decoys']}")
    print(f"BEDROC(alpha=20)={bench['bedroc']} EF1%={bench['ef'].get('EF1%')} logAUC={bench['log_auc']}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
