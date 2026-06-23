"""proteinttt_node — improve low-confidence ESMFold folds via test-time training (§2.3).

ProteinTTT (test-time training) adapts a pretrained protein language model to a
SINGLE target sequence at inference: a few gradient steps of the ESM masked-LM
objective on that sequence (and any homologs) before predicting structure. The
payoff is largest exactly where VTA-Agent is weakest — proteins under-represented
in the training set, i.e. orphan/viral segments (TiLV segments, novel ORFs) whose
ESMFold mean pLDDT lands below the confidence line.

Where it sits:  structure → **proteinttt** → pockets
We only spend the (GPU-heavy) refinement on structures that actually need it:
ESMFold folds whose mean pLDDT is below `_LOW_CONFIDENCE`. Experimental chains
(trusted ground truth, pLDDT=None) and AlphaFold-DB / Boltz-2 models are left
untouched — TTT is specifically an ESM-based-prediction improvement, and a curated
AlphaFold model shouldn't be perturbed by it.

ANNOTATION/REFINEMENT contract (same discipline as the rest of the pipeline):
on success it replaces the low-confidence fold in-place with the improved one
(method → "esmfold+ttt", new mean_plddt, original kept as `refined_from_plddt`);
on absence/failure it leaves the structure exactly as ESMFold produced it. The
graph never crashes for a missing GPU or package.

Auto-detects the `proteinttt` package; degrades to a labelled skip when absent
(true on this arm64/no-GPU box — real refinement runs on a GPU host).

NOTE ON API/CITATION: ProteinTTT post-dates the training cutoff. The
`_refine_structure()` stub reflects the most likely API from the method
description (sequence + initial fold → refined fold + per-residue pLDDT); verify
against the released package and update the call + citation if the signature
differs, exactly as flagged for Boltzina in `vta/nodes/boltzina.py`.
"""
from __future__ import annotations

import os

from vta.state import VTAState

# Same line the structure node uses to print "(LOW CONFIDENCE)" — only refine below it.
_LOW_CONFIDENCE = 70.0

_OUT_DIR = "structures"


def _proteinttt_available() -> bool:
    try:
        import proteinttt  # noqa: F401
        return True
    except Exception:
        return False


def _refine_structure(sequence: str, initial_pdb: str, out_path: str) -> tuple[str, float] | None:
    """Run test-time training to refine one fold; return (pdb_path, mean_plddt) or None.

    Verify the exact API against the released package. The method description is:
    adapt the ESM model to `sequence` at inference, re-predict the structure
    (seeded by the initial ESMFold fold), and write an improved PDB whose B-factor
    column again holds per-residue pLDDT (so `structure.mean_plddt` applies).
    """
    try:
        from proteinttt import refine
        result = refine(sequence=sequence, initial_pdb=initial_pdb, out_pdb=out_path)
        # result may be the pLDDT float, a (path, plddt) pair, or a dict.
        from vta.nodes.structure import mean_plddt
        if isinstance(result, tuple) and len(result) == 2:
            path, plddt = result
            return str(path), float(plddt)
        if isinstance(result, dict):
            path = result.get("pdb_path", out_path)
            plddt = result.get("mean_plddt")
            return str(path), float(plddt) if plddt is not None else mean_plddt(open(path).read())
        if not os.path.exists(out_path):
            return None
        return out_path, mean_plddt(open(out_path).read())
    except Exception:
        return None


def proteinttt_node(state: VTAState) -> VTAState:
    """Refine low-confidence ESMFold structures in place via test-time training."""
    structures = state.get("structures") or {}
    if not structures:
        return state

    # Which folds are eligible: ESMFold output below the confidence line, with a file.
    candidates = [
        name for name, rec in structures.items()
        if rec.get("method") == "esmfold"
        and rec.get("pdb_path")
        and rec.get("mean_plddt") is not None
        and rec["mean_plddt"] < _LOW_CONFIDENCE
    ]

    if not _proteinttt_available():
        if candidates:
            state["audit_trail"].append(
                f"[skip] ProteinTTT: package not installed (pip install proteinttt); "
                f"{len(candidates)} low-confidence fold(s) left as ESMFold produced them")
            state["versions"]["proteinttt"] = "skipped"
        return state

    if not candidates:
        state["audit_trail"].append(
            "ProteinTTT: no low-confidence ESMFold folds to refine (nothing below "
            f"pLDDT {_LOW_CONFIDENCE:.0f})")
        state["versions"]["proteinttt"] = "no-op (all folds above confidence line)"
        return state

    os.makedirs(_OUT_DIR, exist_ok=True)
    refined = 0
    for name in candidates:
        rec = structures[name]
        seq = (state.get("extracted_proteins") or {}).get(name, {}).get("sequence", "")
        out_path = os.path.join(_OUT_DIR, f"{state['run_id']}_{name}_ttt.pdb")
        result = _refine_structure(seq, rec["pdb_path"], out_path)
        if not result:
            state["audit_trail"].append(
                f"ProteinTTT[{name}]: refinement failed; keeping ESMFold fold "
                f"(pLDDT {rec['mean_plddt']})")
            continue
        new_path, new_plddt = result
        old_plddt = rec["mean_plddt"]
        rec.update(
            pdb_path=new_path, mean_plddt=round(new_plddt, 1), method="esmfold+ttt",
            source="ESMFold + ProteinTTT (test-time training)",
            refined_from_plddt=old_plddt,
        )
        refined += 1
        gain = round(new_plddt - old_plddt, 1)
        flag = " (still LOW CONFIDENCE)" if new_plddt < _LOW_CONFIDENCE else ""
        state["audit_trail"].append(
            f"ProteinTTT[{name}]: refined ESMFold fold, pLDDT {old_plddt} → "
            f"{round(new_plddt, 1)} ({'+' if gain >= 0 else ''}{gain}){flag}")

    state["versions"]["proteinttt"] = (
        "ProteinTTT (test-time training)" if refined
        else "skipped (refinement unavailable for all candidates)")
    return state
