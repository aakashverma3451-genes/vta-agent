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
    candidate_targets: Optional[List[Dict]]
    target_prioritization: Optional[List[Dict]]

    # --- Module 2: structure folding --------------------------------------
    structures: Optional[Dict[str, Dict]]      # {name: {pdb_path, mean_plddt}}

    # --- Module 2b: pockets (Phase 2; mocked now) -------------------------
    pockets: Optional[Dict[str, List[Dict]]]   # {protein: [{id, center, druggability, ...}]}
    # Per-residue conservation (JSD) keyed by resseq, written by conservation_node
    # alongside the pocket aggregate. conservation_contacts_node (SPEC #4) reads it to
    # weight each ligand by the conservation of the residues ITS pose contacts.
    residue_conservation: Optional[Dict[str, Dict]]  # {protein: {resseq: jsd}}
    docking_species: Optional[Dict[str, Dict]]        # {ligand: resolved docking species}

    # --- Module 3: docking (mocked in Phase 1) ----------------------------
    # docking records may gain cnn_score/cnn_affinity from the DL-rescore seam
    # (rescore_node, Phase 3b) when gnina is present — annotation-only, ranking
    # is unchanged until the term is calibrated against the validation gate.
    docking_results: Optional[List[Dict]]      # [{ligand, pocket, dG, rmsd, le, ...}]
    lead_candidates: Optional[List[Dict]]      # top-N ranked
    chemistry_annotations: Optional[Dict[str, Dict]]  # {ligand: {chemistry, active_species}}
    counter_target_panel: Optional[List[Dict]]
    resistance_mutants: Optional[List[Dict]]
    prediction_registry: Optional[Dict[str, Any]]

    # --- Phase 4: MD validation (opt-in, include_md=True) ----------------
    # Pocket records gain consensus:bool + detectors:[str] from P2Rank (§3.1).
    md_candidates: Optional[List[Dict]]         # top-N leads selected for MD
    md_results: Optional[Dict[str, Dict]]       # {ligand: {trajectory, status, ...}}
    md_analysis: Optional[Dict[str, Dict]]      # {ligand: {rmsd, contacts, mmgbsa, verdict}}
    md_validated_leads: Optional[List[Dict]]    # md_rerank output (md_score, md_badge)

    # --- Phase 5: FEP / ABFE validation (opt-in, include_fep=True) --------
    fep_results: Optional[Dict[str, Dict]]       # {ligand: {delta_g, error, status, ...}}
    fep_validated_leads: Optional[List[Dict]]    # top leads annotated with fep_* fields

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
