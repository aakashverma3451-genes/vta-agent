"""Stratify the Mpro dataset by binding mode (covalent vs non-covalent).

AutoDock Vina scores non-covalent binding only, so a covalent inhibitor's affinity is not
meaningful in a non-covalent dock. The headline powered benchmark therefore uses
**non-covalent actives only**; covalent actives are kept as a separate, clearly-labelled
sub-list and never pooled — the same binding-mode discipline that separated HCV NI from NNI.

Inactives are controls regardless of warhead, but each is annotated with its binding mode
for transparency. The non-covalent active list written here is the benchmark's active input.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from vta.chem.warheads import classify_binding_mode

DATASET = Path("vta/data/mpro/mpro_dataset.json")
OUT = Path("vta/data/mpro/mpro_stratified.json")


def _classify(records: list[dict]) -> list[dict]:
    out = []
    for r in records:
        mode = classify_binding_mode(r.get("smiles"), r.get("covalent_annotation") or None)
        out.append({**r, "binding_mode": mode["mode"], "binding_mode_basis": mode["basis"],
                    "warheads": mode["warheads"]})
    return out


def stratify(dataset: dict) -> dict:
    actives = _classify(dataset.get("actives", []))
    inactives = _classify(dataset.get("inactives", []))

    def _split(records: list[dict], mode: str) -> list[dict]:
        return [r for r in records if r["binding_mode"] == mode]

    non_covalent_actives = _split(actives, "non_covalent")
    covalent_actives = _split(actives, "covalent")
    ambiguous_actives = _split(actives, "ambiguous")

    def _mode_counts(records: list[dict]) -> dict:
        return dict(Counter(r["binding_mode"] for r in records))

    def _source_counts(records: list[dict]) -> dict:
        return dict(Counter(r["source"] for r in records))

    return {
        "target": dataset.get("target"),
        "binding_mode_rule": (
            "Covalent = explicit Moonshot/ChEMBL covalent flag OR Cys145 electrophilic "
            "warhead SMARTS match. Headline benchmark uses non_covalent actives only; "
            "covalent actives are a separate sub-analysis (not pooled)."
        ),
        "active_mode_counts": _mode_counts(actives),
        "inactive_mode_counts": _mode_counts(inactives),
        "non_covalent_active_sources": _source_counts(non_covalent_actives),
        "inactive_sources": _source_counts(inactives),
        "n_non_covalent_actives": len(non_covalent_actives),
        "n_covalent_actives": len(covalent_actives),
        "n_ambiguous_actives": len(ambiguous_actives),
        "n_inactives": len(inactives),
        "headline_benchmark_inputs": {
            "actives": "non_covalent_actives",
            "controls": "inactives (experimentally measured)",
            "note": "covalent_actives held out; require covalent docking, not Vina.",
        },
        "non_covalent_actives": non_covalent_actives,
        "covalent_actives": covalent_actives,
        "ambiguous_actives": ambiguous_actives,
        "inactives": inactives,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", default=str(DATASET))
    args = parser.parse_args()
    dataset = json.loads(Path(args.dataset).read_text())
    strat = stratify(dataset)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(strat, indent=2))
    print(f"active mode counts: {strat['active_mode_counts']}")
    print(f"non-covalent actives: {strat['n_non_covalent_actives']} "
          f"(sources {strat['non_covalent_active_sources']})")
    print(f"covalent actives (held out): {strat['n_covalent_actives']}")
    print(f"inactives (controls): {strat['n_inactives']} (sources {strat['inactive_sources']})")
    print(OUT)


if __name__ == "__main__":
    main()
