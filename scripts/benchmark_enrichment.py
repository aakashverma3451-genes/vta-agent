"""Retrospective enrichment benchmark for the VTA-Agent screen (SPEC #5).

SUPPLEMENTS `scripts/validate_controls.py` (the fast 9-ligand pre-commit sanity check)
with the enrichment statistics a paper needs: dock a curated ACTIVE set (nucleotide-
analog RdRp inhibitors) plus DECOYS into the PB1 NTP site, rank by score, and quantify
separation with EF1%/EF5%/EF10%/EF20%, BEDROC (α=20, early recognition), and ROC-AUC.

Score source, in priority order (honest provenance, never fabricated):
  1. a committed REAL-score cache (`vta/data/benchmark_scores.json`, from a real Vina
     run via validate_controls.py) — used here so the benchmark reproduces offline;
  2. (extensible) a live REAL chain over actives_rdrp.smi + fetch_decoys when VINA_BIN +
     network are present — the path to a larger, rigorous benchmark.

HONESTY: the committed run is SMALL N (the gate set) and its decoys are mechanism-
distinct antivirals, NOT DUD-E property-matched. DUD-E has documented analog bias that
inflates scores; LIT-PCBA is the less-biased alternative (a 2025 audit found leakage even
there). Treat these numbers as a demonstration/baseline, not a published result — see
ARCHITECTURE.md §8 (Data leakage & evaluation honesty).

Run from the vta-agent dir:
    python scripts/benchmark_enrichment.py
"""
from __future__ import annotations

import json
import os

from vta.eval.metrics import enrichment_report

_SCORE_CACHE = os.path.join("vta", "data", "benchmark_scores.json")
_OUT = os.path.join("outputs", "benchmark_enrichment.json")

_CAVEAT = (
    "CAVEAT (read before quoting): small N (gate set); decoys are mechanism-distinct "
    "antivirals, NOT DUD-E property-matched. DUD-E has analog bias that inflates scores "
    "(LIT-PCBA is the less-biased alternative; a 2025 audit found leakage even there). "
    "This is a demonstration/baseline — see ARCHITECTURE.md §8."
)


def load_scored_entries(cache_path: str = _SCORE_CACHE) -> tuple[list[dict], str]:
    """Load committed real-score entries; return (entries, provenance)."""
    with open(cache_path) as fh:
        data = json.load(fh)
    return data.get("entries", []), data.get("_provenance", "(no provenance recorded)")


def main() -> None:
    if not os.path.exists(_SCORE_CACHE):
        print(f"[skip] no score cache at {_SCORE_CACHE}; run scripts/validate_controls.py "
              "with VINA_BIN to populate real scores first.")
        return

    entries, provenance = load_scored_entries()
    rep = enrichment_report(entries)

    print("\n=== VTA-Agent retrospective enrichment benchmark (SPEC #5) ===")
    print(f"provenance: {provenance}\n")
    print(f"{'rank':<5}{'compound':<22}{'score':>9}  label")
    print("-" * 48)
    for row in rep["ranking"]:
        tag = "★ ACTIVE" if row["label"] == "active" else "decoy"
        print(f"{row['rank']:<5}{row['name']:<22}{row['score']:>9}  {tag}")
    print("-" * 48)
    print(f"N = {rep['n']}  ({rep['n_actives']} actives / {rep['n_decoys']} decoys)")
    ef = "  ".join(f"{k}={v}" for k, v in rep["ef"].items())
    print(f"{ef}")
    print(f"BEDROC(α=20) = {rep['bedroc']}    ROC-AUC = {rep['roc_auc']}")
    print(f"\n{_CAVEAT}")

    os.makedirs("outputs", exist_ok=True)
    with open(_OUT, "w") as fh:
        json.dump({**rep, "_provenance": provenance, "_caveat": _CAVEAT}, fh, indent=2)
    print(f"\nwrote {_OUT}")


if __name__ == "__main__":
    main()
