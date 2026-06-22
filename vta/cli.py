"""vta.cli — the autonomous entry point.

One command takes a viral genome FASTA and runs the full drug-design pipeline by
itself: classify → route → fold → find the catalytic site → dock a real ligand
library → rank → emit an HTML report. No human setup — tools self-configure
(vta.toolconfig) and TaxonAgent is located automatically. The confidence router is
the agent's autonomy: it decides proceed / flag / defer without asking.

    python -m vta run path/to/genome.fasta
    vta run path/to/genome.fasta            # after `pip install -e .`
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _ensure_taxonagent() -> None:
    """Make THIS project's TaxonAgent authoritative without the human setting PYTHONPATH.

    The venv may have a *different* `taxonagent` installed (an older checkout lacking
    the classify/validate_contract re-exports). We prepend our sibling source so it
    wins, and evict any already-imported wrong one — otherwise the autonomous run
    silently picks up the stale package and fails at the classify node.
    """
    for parents_up in (2, 3):
        cand = Path(__file__).resolve().parents[parents_up] / "taxonagent" / "src"
        if (cand / "taxonagent" / "contract.py").exists():   # our version (has classify)
            sys.path.insert(0, str(cand))
            mod = sys.modules.get("taxonagent")
            if mod and not str(getattr(mod, "__file__", "")).startswith(str(cand)):
                for k in [m for m in list(sys.modules)
                          if m == "taxonagent" or m.startswith("taxonagent.")]:
                    del sys.modules[k]
            return


def _summary(state: dict) -> None:
    tr = state.get("taxon_result") or {}
    print("\n" + "=" * 60)
    print(f"  {tr.get('species', '?')}  ({tr.get('genus', '?')})")
    print(f"  confidence {state.get('classification_confidence')}%  "
          f"-> {str(state.get('route', '?')).upper()}")
    leads = state.get("lead_candidates") or []
    if leads:
        print(f"\n  Top leads ({len(leads)}):")
        for i, r in enumerate(leads[:5], 1):
            star = " *" if r.get("positive_control") else "  "
            print(f"   {i}.{star} {r.get('ligand'):<18} "
                  f"dG={r.get('dG')}  LE={r.get('le')}  score={r.get('score')}")
    else:
        print("\n  No leads — deferred to a human expert.")
    rep = next((l for l in state.get("audit_trail", []) if l.startswith("Report:")), None)
    if rep:
        print(f"\n  {rep}")
    print("=" * 60)


def run(fasta: str, run_id: str | None = None) -> dict:
    """Run the full pipeline on a genome FASTA and print a summary. Returns state."""
    _ensure_taxonagent()
    from vta.graph import build_app
    from vta.state import new_state
    from vta.toolconfig import find_tool

    rid = run_id or Path(fasta).stem
    print(f"[vta] genome : {fasta}")
    print(f"[vta] fpocket: {find_tool('fpocket', 'FPOCKET_BIN') or 'MOCK (not found)'}")
    print(f"[vta] vina   : {find_tool('vina', 'VINA_BIN') or 'MOCK (not found)'}")
    print("[vta] running autonomous pipeline (real docking can take a few minutes)…")

    final = build_app().invoke(new_state(fasta, rid))
    _summary(final)
    return final


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="vta", description="Autonomous viral target assessment (genome -> drug leads)")
    sub = p.add_subparsers(dest="cmd")
    rp = sub.add_parser("run", help="run the pipeline on a genome FASTA")
    rp.add_argument("fasta", help="path to a viral genome FASTA")
    rp.add_argument("--run-id", default=None, help="label for this run (default: filename)")

    args = p.parse_args(argv)
    if args.cmd == "run":
        run(args.fasta, args.run_id)
        return 0
    p.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
