"""triage_router_node (Phase R / R2) — the central lesson: route each target to the right method.

HemaGuide's measured "critical integrating element" was autonomous routing to case-specific
reasoning modes (plain LLM 0% → full agent 86.7%). The VTA analogue routes each *target* to
one of {full_dock | annotate_only | defer | refuse} from its dossier (R1), instead of applying
rigid Vina to everything. This is where the honest value concentrates: it stops the pipeline
producing a spuriously-precise docking ranking for targets where docking is out-of-domain.

Design invariant — this changes WHICH method runs, never the docking weights/scoring.

Routing is deterministic (rules below; an LLM may later only *summarise* the rationale, never
override it) and **defaults to full_dock**. A downgrade fires only on positive, run-available
evidence, so a drug-like screen against any target still docks:

  refuse       — no usable structure (provenance refused, or binding-site pLDDT below a hard floor)
  defer        — borderline predicted structure (flagged QC / intermediate binding-site pLDDT)
  annotate_only— scoring is out-of-domain: benchmarkability un_benchmarkable, OR a metal ion is
                 in the *detected pocket*, OR the screened ligand class is charged/nucleotide/
                 covalent, OR the R3 playbook prior says rigid Vina lost to the 2D baseline here
  full_dock    — otherwise (experimental/adequate pocket + in-domain ligands + benchmarkable)

Note the router keys on `pocket_descriptors.metal_ions_present` (a metal actually detected in
THIS run's pocket) and on benchmarkability, NOT on the target's abstract `metal_dependence` /
`known_ligand_class` — so a metal-dependent polymerase screened with drug-like compounds and no
metal resolved in its pocket still docks (e.g. the TiLV PB1 demonstration).
"""
from __future__ import annotations

from typing import Any, Dict, List

from vta.state import VTAState

# binding-site pLDDT gates (only applied to PREDICTED structures with a numeric proxy value;
# experimental structures carry None and skip these). Mirrors structure_qc's <70 low-pLDDT flag.
PLDDT_REFUSE_BELOW = 50.0
PLDDT_DEFER_BELOW = 70.0
_OUT_OF_DOMAIN_LIGANDS = {"charged", "nucleotide", "covalent"}


def _decide(dossier: Dict[str, Any], playbook: Dict[str, Any], screened_ligand_class):
    """Return (decision, rationale[]) for one target from its dossier + R3 prior + ligand hint."""
    rationale: List[str] = []
    sq = dossier.get("structure_quality") or {}
    pocket = dossier.get("pocket_descriptors") or {}
    bench = (dossier.get("benchmarkability") or {}).get("status")
    provenance = sq.get("provenance")
    plddt = sq.get("binding_site_plddt")
    predicted = provenance == "predicted"

    # 1) refuse — no usable structure
    if provenance == "refused" or provenance == "none" or not sq.get("method"):
        return "refuse", ["no usable structure (provenance=%s)" % provenance]
    if predicted and isinstance(plddt, (int, float)) and plddt < PLDDT_REFUSE_BELOW:
        return "refuse", [f"binding-site pLDDT {plddt} < {PLDDT_REFUSE_BELOW} (structure unusable)"]

    # 2) defer — borderline predicted structure
    if predicted and isinstance(plddt, (int, float)) and plddt < PLDDT_DEFER_BELOW:
        return "defer", [f"predicted structure, binding-site pLDDT {plddt} in "
                         f"[{PLDDT_REFUSE_BELOW}, {PLDDT_DEFER_BELOW}) — borderline"]
    if predicted and (sq.get("qc_status") == "flagged"):
        return "defer", ["predicted structure flagged by QC — borderline docking-appropriateness"]

    # 3) annotate_only — scoring out-of-domain (positive, run-available evidence only)
    if bench == "un_benchmarkable":
        rationale.append("benchmarkability=un_benchmarkable (matched-decoy validation infeasible)")
    if pocket.get("metal_ions_present"):
        rationale.append("metal ion resolved in the docked pocket (scoring out-of-domain)")
    if screened_ligand_class in _OUT_OF_DOMAIN_LIGANDS:
        rationale.append(f"screened ligand class '{screened_ligand_class}' is out-of-domain for Vina")
    if (playbook or {}).get("status") == "lost_to_2d_baseline":
        rationale.append("R3 playbook: rigid Vina lost to the 2D baseline for this pocket class")
    if rationale:
        return "annotate_only", rationale

    # 4) full_dock — default
    return "full_dock", ["experimental/adequate structure, in-domain ligands, no downgrade trigger"]


def triage_router_node(state: VTAState) -> VTAState:
    """Route each dossiered target to full_dock | annotate_only | defer | refuse."""
    dossier = state.get("target_dossier") or {}
    playbook = state.get("playbook_prior") or {}          # R3 stub: {} → "no precedent"
    # optional per-protein hint about the class of compounds actually being screened
    ligand_classes = state.get("screened_ligand_class") or {}

    decisions: Dict[str, Dict[str, Any]] = {}
    counts: Dict[str, int] = {}
    for protein, dsr in dossier.items():
        decision, rationale = _decide(dsr, playbook.get(protein), ligand_classes.get(protein))
        decisions[protein] = {"decision": decision, "rationale": rationale}
        counts[decision] = counts.get(decision, 0) + 1

    state["triage_decision"] = decisions
    # keep the R3 stub explicit so downstream/audit sees grounding status
    if not state.get("playbook_prior"):
        state["playbook_prior"] = {p: {"status": "no_precedent",
                                       "note": "R3 PlaybookMemory not populated"} for p in dossier}
    summary = ", ".join(f"{k}={v}" for k, v in sorted(counts.items())) or "none"
    state["audit_trail"].append(f"Triage: routed {len(decisions)} target(s) [{summary}]")
    state["versions"]["triage"] = "triage-router v1 (R2); R3 playbook stubbed"
    return state
