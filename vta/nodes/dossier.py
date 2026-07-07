"""dossier_node (Phase R / R1) — assemble a structured target dossier before any routing.

HemaGuide converts messy input into a section-aware structured case *before* reasoning; the
measured lesson is that making the input legible is what lets the router pick the right method
per case. The VTA analogue is a typed **target dossier**, one per protein, assembled from the
structure/pocket/species nodes plus a small curated target-knowledge table — read by the
TriageRouter (R2) to decide full_dock | annotate_only | defer | refuse.

Deterministic and honest: every field comes from state already present or from a committed
knowledge table (domain facts, mirroring the EXPERIMENTAL_PDB tables — not fabricated). A field
we cannot ground is `null` with a reason, never invented. In particular, per-binding-site pLDDT
is NOT wired in this codebase, so `binding_site_plddt` carries the mean-pLDDT proxy with an
explicit note (the Scardino/Cavasotto 2023 caveat: global pLDDT is necessary-not-sufficient for
docking-appropriateness — recorded so the router and report never over-trust it).
"""
from __future__ import annotations

from typing import Any, Dict, Optional

from vta.report_envelope import DEFAULT_BENCHMARK_PATH, _load_json, _match_target
from vta.state import VTAState

# Curated target biology (domain knowledge, not fabrication). Keyword-matched against protein
# names. metal_dependence and typical_ligand_class are the target's *biology*; they are recorded
# for context but do NOT by themselves trigger a docking downgrade — the router keys on
# run-available evidence (metal actually in the detected pocket, benchmarkability, pocket-site
# confidence, the screened ligand class), so a drug-like screen against a polymerase still docks.
_TARGET_KNOWLEDGE = (
    (("mpro", "3clpro", "main protease"),
     {"target_class": "protease", "catalytic_mechanism": "cysteine protease (His41/Cys145 dyad)",
      "metal_dependence": False, "typical_ligand_class": "drug_like"}),
    (("ns5b", "hcv"),
     {"target_class": "polymerase (RdRp)", "catalytic_mechanism": "two-metal nucleotidyl transfer",
      "metal_dependence": True, "typical_ligand_class": "nucleotide"}),
    (("pb1", "tilv", "rdrp", "polymerase"),
     {"target_class": "polymerase (RdRp)", "catalytic_mechanism": "two-metal nucleotidyl transfer",
      "metal_dependence": True, "typical_ligand_class": "nucleotide"}),
)

_PREDICTED = {"esmfold", "alphafold", "boltz2", "boltz-2", "proteinttt"}


def _knowledge(protein: str) -> Dict[str, Any]:
    p = (protein or "").lower()
    for keywords, facts in _TARGET_KNOWLEDGE:
        if any(k in p for k in keywords):
            return dict(facts)
    return {"target_class": None, "catalytic_mechanism": None,
            "metal_dependence": None, "typical_ligand_class": None}


def _benchmarkability(protein: str, benchmark: Optional[dict]) -> Dict[str, Any]:
    """Grade → {powered | underpowered | un_benchmarkable | unknown} from the frozen benchmark."""
    matched = _match_target(protein)
    if not matched:
        return {"status": "unknown", "reason": "no frozen benchmark covers this target"}
    frozen_name = matched[0]
    for t in ((benchmark or {}).get("targets") or []):
        if t.get("target") == frozen_name:
            grade = (t.get("grade") or "").lower()
            if "not benchmarkable" in grade or "un_benchmarkable" in grade:
                status = "un_benchmarkable"
            elif "underpowered" in grade or "under-powered" in grade:
                status = "underpowered"
            elif "powered" in grade:
                status = "powered"
            else:
                status = "unknown"
            return {"status": status, "reason": t.get("grade"), "benchmark_target": frozen_name}
    return {"status": "unknown", "reason": "target not present in frozen benchmark"}


def _pocket_descriptors(pockets: Optional[list]) -> Dict[str, Any]:
    if not pockets:
        return {"available": False, "metal_ions_present": False}
    top = pockets[0]
    metals = top.get("metal_ions") or top.get("metal_ions_present")
    return {
        "available": True,
        "source": top.get("source") or top.get("detector") or "detected",
        "volume": top.get("volume"),
        "druggability": top.get("druggability"),
        "center": top.get("center"),
        # metal actually detected in THIS run's pocket (not the target's abstract biology)
        "metal_ions_present": bool(metals),
        "metal_ions": metals if isinstance(metals, list) else None,
    }


def _structure_quality(record: Dict[str, Any]) -> Dict[str, Any]:
    method = (record.get("method") or "none").lower()
    provenance = ("experimental" if method == "experimental"
                  else "predicted" if method in _PREDICTED
                  else "refused" if "refuse" in method
                  else method)
    qc = record.get("structure_qc") or {}
    return {
        "provenance": provenance,
        "method": method,
        "resolution": record.get("resolution"),
        # per-binding-site pLDDT is not wired; carry mean pLDDT as a labelled proxy
        "binding_site_plddt": record.get("mean_plddt"),
        "binding_site_plddt_note": ("proxy: mean pLDDT — per-binding-site pLDDT not wired "
                                    "(global pLDDT is necessary-not-sufficient; Scardino 2023)"),
        "apo_holo": record.get("apo_holo") or "unknown",
        "qc_status": qc.get("status", "unknown"),
        "qc_issues": qc.get("issues") or qc.get("issue") or [],
    }


def dossier_node(state: VTAState) -> VTAState:
    """Build one target dossier per protein with a structure, before routing."""
    structures = state.get("structures") or {}
    pockets = state.get("pockets") or {}
    benchmark = _load_json(DEFAULT_BENCHMARK_PATH)

    dossier: Dict[str, Dict[str, Any]] = {}
    for name, record in structures.items():
        know = _knowledge(name)
        dossier[name] = {
            "protein": name,
            "structure_quality": _structure_quality(record),
            "pocket_descriptors": _pocket_descriptors(pockets.get(name)),
            "target_class": know["target_class"],
            "catalytic_mechanism": know["catalytic_mechanism"],
            "metal_dependence": know["metal_dependence"],
            "known_ligand_class": know["typical_ligand_class"],
            "benchmarkability": _benchmarkability(name, benchmark),
        }

    state["target_dossier"] = dossier
    state["audit_trail"].append(
        f"Dossier: assembled {len(dossier)} target dossier(s) "
        f"[{', '.join(sorted(dossier)) or 'none'}]")
    state["versions"]["dossier"] = "target-dossier v1 (R1)"
    return state
