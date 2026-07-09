"""vta.cli — the autonomous entry point.

One command takes a viral genome FASTA and runs the full drug-design pipeline by
itself: classify → route → fold → find the catalytic site → dock a real ligand
library → rank → emit an HTML report. Tools self-configure through `vta.toolconfig`;
TaxonAgent is imported as a normal installed package. The confidence router is the
agent's autonomy: it decides proceed / flag / defer without asking.

    python -m vta run path/to/genome.fasta
    vta run path/to/genome.fasta            # after `pip install -e .`
"""
from __future__ import annotations

import argparse
from pathlib import Path


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
    elif state.get("route") == "defer":
        print("\n  No leads — confidence too low, deferred to a human expert.")
    else:
        print("\n  No leads — proceeded but docking produced no ranked candidates "
              "(see audit trail).")
    rep = next((l for l in state.get("audit_trail", []) if l.startswith("Report:")), None)
    if rep:
        print(f"\n  {rep}")
    print("=" * 60)


def run(fasta: str, run_id: str | None = None, *, include_md: bool = False,
        include_fep: bool = False) -> dict:
    """Run the full pipeline on a genome FASTA and print a summary. Returns state."""
    from vta.graph import build_app
    from vta.state import new_state
    from vta.toolconfig import find_tool

    rid = run_id or Path(fasta).stem
    print(f"[vta] genome : {fasta}")
    print(f"[vta] fpocket: {find_tool('fpocket', 'FPOCKET_BIN') or 'MOCK (not found)'}")
    print(f"[vta] vina   : {find_tool('vina', 'VINA_BIN') or 'MOCK (not found)'}")
    print("[vta] running autonomous pipeline (real docking can take a few minutes)…")

    final = build_app(include_md=include_md, include_fep=include_fep).invoke(
        new_state(fasta, rid))
    _summary(final)
    return final


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="vta", description="Autonomous viral target assessment (genome -> drug leads)")
    sub = p.add_subparsers(dest="cmd")
    rp = sub.add_parser("run", help="run the pipeline on a genome FASTA")
    rp.add_argument("fasta", help="path to a viral genome FASTA")
    rp.add_argument("--run-id", default=None, help="label for this run (default: filename)")
    rp.add_argument("--include-md", action="store_true",
                    help="run the opt-in MD validation phase")
    rp.add_argument("--include-fep", action="store_true",
                    help="run opt-in FEP/ABFE after MD validation")

    args = p.parse_args(argv)
    if args.cmd == "run":
        run(args.fasta, args.run_id, include_md=args.include_md,
            include_fep=args.include_fep)
        return 0
    p.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
