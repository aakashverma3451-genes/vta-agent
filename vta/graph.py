"""vta.graph — compile the VTA-Agent LangGraph application.

Three build modes controlled by `build_app(include_md=False, include_fep=False)`:

  Phase A (default, fast):
    classify → route → structure → proteinttt → pockets → conservation → species_resolution → dock → rescore → rank → chemistry → admet → report

  Phase B (opt-in, slow — requires OpenMM + MDAnalysis):
    …same up to admet… → md_select → md_simulate → md_analyze → md_rerank → report

  Phase C (opt-in, very slow — requires an external ABFE runner):
    …same up to md_rerank… → fep → report

Both proceed paths and the defer path always end at report_node so every run
emits a self-contained HTML report.
"""
from __future__ import annotations

from langgraph.graph import END, StateGraph

from vta.nodes.classify import classify_node
from vta.nodes.conservation import conservation_node
from vta.nodes.conservation_contacts import conservation_contacts_node
from vta.nodes.consensus import consensus_node
from vta.nodes.docking import docking_node
from vta.nodes.dossier import dossier_node
from vta.nodes.pockets import pockets_node
from vta.nodes.proteinttt import proteinttt_node
from vta.nodes.admet import admet_node
from vta.nodes.chemistry import chemistry_node
from vta.nodes.md_select import md_select_node
from vta.nodes.md_simulate import md_simulate_node
from vta.nodes.md_analyze import md_analyze_node
from vta.nodes.md_rerank import md_rerank_node
from vta.nodes.fep import fep_node
from vta.nodes.boltzina import boltzina_node
from vta.nodes.rank import rank_node
from vta.nodes.rescore import rescore_node
from vta.nodes.resistance import resistance_node
from vta.nodes.report import report_node
from vta.nodes.router import (
    DEFER_NODE,
    STRUCTURE_NODE,
    route_by_confidence,
    route_edge,
)
from vta.nodes.species_resolution import species_resolution_node
from vta.nodes.selectivity import selectivity_node
from vta.nodes.stubs import defer_node
from vta.nodes.structure import structure_node
from vta.nodes.structure_qc import structure_qc_node
from vta.nodes.target_prioritization import target_prioritization_node
from vta.nodes.triage import triage_router_node
from vta.state import VTAState


def build_app(include_md: bool = False, include_fep: bool = False):
    """Compile the VTA-Agent graph.

    Args:
        include_md: If True, append the MD validation phase after admet
                    (md_select → md_simulate → md_analyze → md_rerank).
                    MD is slow (GPU-hours for production) so it's opt-in;
                    the default fast path skips it and goes straight to report.
        include_fep: If True, append FEP/ABFE validation after MD re-ranking.
                     This implies include_md because FEP validates top MD leads.
    """
    include_md = include_md or include_fep
    g = StateGraph(VTAState)

    # ── shared nodes (both modes) ────────────────────────────────────────
    g.add_node("classify", classify_node)
    g.add_node("router", route_by_confidence)
    g.add_node(STRUCTURE_NODE, structure_node)
    g.add_node("structure_qc", structure_qc_node)
    g.add_node("target_prioritization", target_prioritization_node)
    g.add_node("proteinttt", proteinttt_node)  # refine low-pLDDT ESMFold folds (§2.3)
    g.add_node("pockets", pockets_node)
    g.add_node("conservation", conservation_node)  # real per-pocket JSD conservation (§SPEC#1)
    g.add_node("species_resolution", species_resolution_node)
    g.add_node("dossier", dossier_node)  # R1: structured target dossier before routing
    g.add_node("triage", triage_router_node)  # R2: route full_dock|annotate_only|defer|refuse
    g.add_node("dock", docking_node)
    g.add_node("conservation_contacts", conservation_contacts_node)  # ligand-weighted JSD (§SPEC#4)
    g.add_node("rescore", rescore_node)    # GNINA CNN re-score (annotation-only)
    g.add_node("boltzina", boltzina_node)  # Boltzina DL affinity seam (annotation-only)
    g.add_node("consensus", consensus_node)
    g.add_node("rank", rank_node)
    g.add_node("chemistry", chemistry_node)
    g.add_node("selectivity", selectivity_node)
    g.add_node("resistance", resistance_node)
    g.add_node("admet", admet_node)
    g.add_node(DEFER_NODE, defer_node)
    g.add_node("report", report_node)

    g.set_entry_point("classify")
    g.add_edge("classify", "router")
    g.add_conditional_edges("router", route_edge, {
        STRUCTURE_NODE: STRUCTURE_NODE,
        DEFER_NODE: DEFER_NODE,
    })
    g.add_edge(STRUCTURE_NODE, "structure_qc")
    g.add_edge("structure_qc", "target_prioritization")
    g.add_edge("target_prioritization", "proteinttt")  # TTT refinement (skips w/o package)
    g.add_edge("proteinttt", "pockets")
    g.add_edge("pockets", "conservation")   # overwrite 0.5 placeholder w/ real JSD (skips w/o MSA)
    g.add_edge("conservation", "species_resolution")
    g.add_edge("species_resolution", "dossier")  # R1 dossier before docking
    g.add_edge("dossier", "triage")              # R2 route the target
    g.add_edge("triage", "dock")
    g.add_edge("dock", "conservation_contacts")   # ligand-contact-weight conservation (§SPEC#4)
    g.add_edge("conservation_contacts", "rescore")  # GNINA CNN re-score (skips w/o gnina binary)
    g.add_edge("rescore", "boltzina")  # Boltzina DL affinity (skips w/o package)
    g.add_edge("boltzina", "consensus")
    g.add_edge("consensus", "rank")
    g.add_edge("rank", "chemistry")
    g.add_edge("chemistry", "selectivity")
    g.add_edge("selectivity", "resistance")
    g.add_edge("resistance", "admet")
    g.add_edge(DEFER_NODE, "report")
    g.add_edge("report", END)

    if include_md:
        # ── Phase B: MD validation (opt-in, hours) ───────────────────────
        g.add_node("md_select", md_select_node)
        g.add_node("md_simulate", md_simulate_node)
        g.add_node("md_analyze", md_analyze_node)
        g.add_node("md_rerank", md_rerank_node)
        g.add_edge("admet", "md_select")
        g.add_edge("md_select", "md_simulate")
        g.add_edge("md_simulate", "md_analyze")
        g.add_edge("md_analyze", "md_rerank")
        if include_fep:
            g.add_node("fep", fep_node)
            g.add_edge("md_rerank", "fep")
            g.add_edge("fep", "report")
        else:
            g.add_edge("md_rerank", "report")
    else:
        # ── Phase A: fast path ───────────────────────────────────────────
        g.add_edge("admet", "report")

    return g.compile()
