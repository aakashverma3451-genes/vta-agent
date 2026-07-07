"""annotate_rank_node (Phase R / R5) — the annotate_only executor.

When the TriageRouter (R2) routes a target away from docking, or the VerificationNode (R4)
downgrades a docking claim, the pipeline must still return something useful — but **honestly
labelled** as a ligand-based annotation, never a docking-precise enrichment number.

This node produces, for each `annotate_only` target, a labelled ligand-based ranking:
  - if a known-actives library is available (state["annotation_library"] with positive_control
    flags), rank candidates by 2-D ECFP4 similarity to the actives (the same trivial baseline
    the validation layer uses) — explicitly tagged `ligand_based_2d_similarity_annotation`;
  - otherwise emit a labelled `annotation_only` record stating no validated ranking is possible
    and why.
It never fabricates compounds and never attaches a ΔG / enrichment claim. Metal-dependent /
nucleotide targets carry the Phase-9D scoring caveat.
"""
from __future__ import annotations

from typing import Any, Dict, List

from vta.eval.baselines import baseline_2d_similarity
from vta.state import VTAState

_NUCLEOTIDE_CAVEAT = (
    "Nucleotide/metal target: docking scores the parent species, not the active "
    "(triphosphate/metal-coordinated) form; this is a ligand-based annotation, NOT a "
    "validated structure-based enrichment (Phase 9D).")


def _annotate_target(protein: str, dossier: Dict[str, Any], library: List[dict]) -> Dict[str, Any]:
    caveats: List[str] = []
    dsr = dossier or {}
    if dsr.get("metal_dependence") or (dsr.get("known_ligand_class") == "nucleotide"):
        caveats.append(_NUCLEOTIDE_CAVEAT)

    actives = [l for l in library if l.get("positive_control")]
    if library and actives:
        smiles = [l.get("smiles") or "" for l in library]
        labels = [1 if l.get("positive_control") else 0 for l in library]
        sim = baseline_2d_similarity(smiles, labels)
        ranked = sorted(
            ({"ligand": l.get("name"), "smiles": l.get("smiles"),
              "similarity_to_actives": round(float(s), 4),
              "positive_control": bool(l.get("positive_control")),
              "ranking_type": "ligand_based_2d_similarity_annotation"}
             for l, s in zip(library, sim)),
            key=lambda r: r["similarity_to_actives"], reverse=True)
        return {"status": "ligand_based_annotation",
                "method": "2D ECFP4 similarity to known actives (leave-one-out)",
                "ranking": ranked, "caveats": caveats,
                "note": "Ligand-based annotation — NOT a docking enrichment claim."}

    return {"status": "annotation_only", "ranking": [], "caveats": caveats,
            "note": ("No validated ranking: docking is out-of-domain for this target and no "
                     "known-actives library was supplied for a 2D-similarity annotation.")}


def annotate_rank_node(state: VTAState) -> VTAState:
    """Produce labelled ligand-based annotations for every annotate_only target."""
    triage = state.get("triage_decision") or {}
    dossier = state.get("target_dossier") or {}
    library = state.get("annotation_library") or []

    annotations: Dict[str, Any] = {}
    for protein, decision in triage.items():
        if (decision or {}).get("decision") == "annotate_only":
            annotations[protein] = _annotate_target(protein, dossier.get(protein), library)

    if annotations:
        state["annotation_rankings"] = annotations
        state["audit_trail"].append(
            f"Annotation ranker: labelled ligand-based annotation for "
            f"{len(annotations)} annotate_only target(s) [{', '.join(sorted(annotations))}]")
        state["versions"]["annotate_rank"] = "ligand-based annotation v1 (R5)"
    return state
