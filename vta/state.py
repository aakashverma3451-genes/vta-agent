"""vta.state — the one shared, typed state object every node reads and writes.

This is the load-bearing contract of the whole system (plan Task 1.1). Defining
it first, before any node logic, forces each module to declare exactly what it
produces and consumes — the discipline that prevents integration hell.

Field naming is deliberately aligned with TaxonAgent's real §4.1.1 contract (the
dict `taxonagent.classify()` emits), so the classify node is a straight copy, not
a translation layer. In particular we carry `classification_basis` alongside
`classification_confidence`, because TaxonAgent's `confidence_pct` is an aa-identity
band, NOT a calibrated probability — the router must know which it's thresholding.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, TypedDict


class VTAState(TypedDict, total=False):
    # --- Input -------------------------------------------------------------
    genome_fasta: str                          # input path or raw FASTA/sequence
    run_id: str                                # unique id for this run (audit)

    # --- Module 1: TaxonAgent output (the §4.1.1 contract) -----------------
    taxon_result: Optional[Dict[str, Any]]     # full TaxonAgent contract dict
    extracted_proteins: Optional[Dict[str, Dict]]   # {name: {sequence, length_aa, plddt, ...}}
    classification_confidence: Optional[float]      # confidence_pct, 0-100
    classification_basis: Optional[str]             # confidence_basis — what that % MEANS

    # --- Routing decision --------------------------------------------------
    route: Optional[str]                       # "proceed" | "flag" | "defer"

    # --- Module 2: structure folding --------------------------------------
    structures: Optional[Dict[str, Dict]]      # {name: {pdb_path, mean_plddt}}

    # --- Module 2b: pockets (Phase 2; mocked now) -------------------------
    pockets: Optional[Dict[str, List[Dict]]]   # {protein: [{id, center, druggability, ...}]}

    # --- Module 3: docking (mocked in Phase 1) ----------------------------
    docking_results: Optional[List[Dict]]      # [{ligand, pocket, dG, rmsd, le, ...}]
    lead_candidates: Optional[List[Dict]]      # top-N ranked

    # --- Cross-cutting: audit + reproducibility ---------------------------
    audit_trail: List[str]                     # every decision, appended in order
    versions: Dict[str, str]                   # tool versions, for reproducibility


def new_state(genome_fasta: str, run_id: str) -> VTAState:
    """Build a fresh state with the cross-cutting accumulators initialised.

    Nodes append to `audit_trail` / `versions`, so they must exist (not be None)
    before the graph runs. Everything else is filled in by the nodes themselves.
    """
    return VTAState(
        genome_fasta=genome_fasta,
        run_id=run_id,
        audit_trail=[],
        versions={},
    )
