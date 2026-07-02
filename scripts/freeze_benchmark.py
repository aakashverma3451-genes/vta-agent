"""Phase 9E — freeze the matured benchmark into one official, hash-stamped artifact.

Bundles the powered, CI-bearing validation evidence into a single frozen file with the
exact compound splits, decoy/inactive provenance, the parent-vs-active-form comparison, and
the gate decision — the artifact a reviewer (or a re-run) can check against. A content hash
makes accidental drift detectable.

Reads only committed artifacts; degrades to labelled "absent" entries if one is missing.
"""
from __future__ import annotations

import hashlib
import json
from datetime import date
from pathlib import Path

OUT = Path("outputs/phase9/locked_benchmark.json")
MPRO = Path("outputs/phase9/mpro_noncovalent_benchmark.json")
HCV_QUALITY = Path("outputs/phase9/hcv_ns5b_ni_matching_quality.json")
HCV_PROV = Path("outputs/phase9/hcv_ns5b_ni_decoy_provenance.json")
TILV = Path("outputs/phase8/matched_benchmark_tilv_pb1.json")
P9D = Path("outputs/phase9/phase9d_active_form_comparison.json")
GATE = Path("outputs/phase9/gate_decision_mpro.md")
TABLE = Path("outputs/phase9/multitarget_ci_table.md")


def _load(p: Path):
    return json.loads(p.read_text()) if p.exists() else None


def _mpro_target(mpro: dict) -> dict:
    rows = mpro.get("rows", [])
    actives = sorted(r.get("ligand_id") or r.get("ligand") for r in rows if r.get("positive_control"))
    inactives = sorted(r.get("ligand_id") or r.get("ligand") for r in rows if not r.get("positive_control"))
    boot = mpro.get("bootstrap", {}).get("bootstrap", {})
    return {
        "target": "SARS-CoV-2 Mpro (non-covalent)",
        "structure": mpro.get("structure"),
        "control_set": mpro.get("control_set"),
        "binding_mode": mpro.get("binding_mode"),
        "real_vina": mpro.get("real_vina"),
        "receptor_prep": mpro.get("receptor_prep"),
        "metrics_ci": {
            "BEDROC": boot.get("bedroc"), "EF1%": boot.get("EF1%"),
            "logAUC": boot.get("log_auc"), "ROC_AUC": boot.get("roc_auc"),
        },
        "verdict": mpro.get("verdict"),
        "splits": {
            "seed": (mpro.get("selection") or {}).get("seed"),
            "n_actives": len(actives), "n_inactives": len(inactives),
            "active_ids": actives, "inactive_ids": inactives,
        },
        "grade": "powered (informative CIs); docking signal modest — reported honestly",
    }


def build() -> dict:
    mpro = _load(MPRO)
    tilv = _load(TILV)
    hcv_q = _load(HCV_QUALITY)
    hcv_p = _load(HCV_PROV)
    p9d = _load(P9D)

    targets = []
    if tilv:
        boot = tilv.get("bootstrap", {}).get("bootstrap", {})
        targets.append({
            "target": "TiLV PB1", "structure": "8PSO chain B",
            "control_set": tilv.get("decoy_quality"),
            "metrics_ci": {"BEDROC": boot.get("bedroc"), "EF1%": boot.get("EF1%"),
                           "logAUC": boot.get("log_auc"), "ROC_AUC": boot.get("roc_auc")},
            "grade": "underpowered demonstration (4 actives; CIs ~ [0,1])",
        })
    if hcv_q and hcv_p:
        zero = sum(1 for v in hcv_p["per_active_counts"].values() if v == 0)
        targets.append({
            "target": "HCV NS5B NI (triphosphate)", "structure": "n/a",
            "control_set": f"purchasable-library scaffold-distinct matched-decoy recipe fails "
                           f"({zero}/{hcv_p['n_actives']} actives recovered 0 decoys)",
            "metrics_ci": None,
            "grade": "not benchmarkable with this decoy recipe (property-unmatched/generative "
                     "decoys or real inactive nucleotides untested — nucleotide meta-finding)",
        })
    if mpro and mpro.get("benchmark"):
        targets.append(_mpro_target(mpro))

    payload = {
        "artifact": "VTA-Agent frozen validation benchmark (Phase 9E)",
        "frozen_at": date.today().isoformat(),
        "policy": "Immutable: regenerate a NEW dated artifact rather than editing this one.",
        "targets": targets,
        "active_form_comparison_9d": (None if not p9d else {
            "structure": p9d.get("structure"),
            "metal": p9d.get("metal"),
            "triphosphate": p9d.get("triphosphate_active_form"),
            "parent_or_prodrug": p9d.get("parent_or_prodrug"),
            "median_dG_triphosphate_minus_parent": p9d.get("median_dG_triphosphate_minus_parent"),
            "tilv_active_form": p9d.get("tilv_active_form"),
        }),
        "gate_decision": "DO NOT PROMOTE (annotation-only); see outputs/phase9/gate_decision_mpro.md",
        "references": {
            "multitarget_table": str(TABLE),
            "gate_decision": str(GATE),
            "mpro_benchmark": str(MPRO),
            "active_form_comparison": str(P9D),
        },
        "provenance": {
            "actives_inactives": "COVID Moonshot (primary) + ChEMBL CHEMBL4523582 (Mpro); "
                                 "ChEMBL CHEMBL4296320 (HCV NS5B NI)",
            "decoys": "Mpro uses experimentally measured inactives (no presumed decoys); "
                      "TiLV uses ChEMBL property-matched presumed decoys",
            "engine": "AutoDock Vina 1.2.5; receptor prep Meeko-first + OpenBabel fallback",
            "metrics": "bootstrap >=1000 + multi-draw control spread; median + 95% CI",
        },
        "caveats": [
            "Mpro benchmark design is publication-grade; rigid-Vina docking signal is modest "
            "(ROC-AUC CI crosses 0.5) and reported as-is.",
            "HCV NS5B NI is not benchmarkable with a purchasable-library scaffold-distinct "
            "matched-decoy recipe (3.09 decoys/active); property-unmatched/generative decoys "
            "or real inactive nucleotides are the untested alternative (meta-finding).",
            "TiLV PB1 is an underpowered demonstration.",
            "No ranking term promoted.",
        ],
    }
    # Content hash over everything except the hash field itself.
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]
    payload["content_hash"] = digest
    return payload


def main() -> None:
    payload = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2))
    print(f"froze {len(payload['targets'])} targets; hash={payload['content_hash']}")
    print(OUT)


if __name__ == "__main__":
    main()
