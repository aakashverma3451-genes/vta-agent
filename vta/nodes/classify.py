"""classify_node — Module 1 seam: run TaxonAgent, validate the contract, record it.

Plan Task 1.2. The node does three things and nothing else:
  1. call the *real* TaxonAgent (imported as a library — we never reimplement it),
  2. validate its output against TaxonAgent's own canonical schema, so a malformed
     or drifted contract fails loudly HERE rather than silently three nodes later,
  3. copy the fields VTA needs into the shared state and log to the audit trail.

#2 (validation) uses `taxonagent.validate_contract` — the single source of truth —
instead of a second hand-rolled schema that could drift from the producer.
"""
from __future__ import annotations

from vta.state import VTAState


def classify_node(state: VTAState) -> VTAState:
    # Real TaxonAgent, re-exported from the package root (cycle-safe, lazy).
    from taxonagent import classify_genome, validate_contract

    result = classify_genome(state["genome_fasta"])
    validate_contract(result)            # fail loud, fail at the seam — not downstream

    state["taxon_result"] = result
    state["extracted_proteins"] = result["extracted_proteins"]
    state["classification_confidence"] = result["confidence_pct"]
    # carry the *meaning* of that %: TaxonAgent's confidence is an aa-identity band,
    # not a calibrated probability — the router must threshold it knowingly.
    state["classification_basis"] = result.get("confidence_basis")
    state["versions"]["taxonagent"] = result["taxonagent_version"]
    state["versions"]["kg"] = result["kg_version"]
    state["audit_trail"].append(
        f"TaxonAgent: {result['genus']} / {result['species']} "
        f"@ {result['confidence_pct']}% (KG {result['kg_version']})"
    )
    return state
