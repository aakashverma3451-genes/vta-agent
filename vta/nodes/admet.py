"""admet_node — autonomous drug-likeness / safety assessment of the leads.

Turns "did it dock?" into "is it a plausible drug?" — the step that makes the
autonomous pipeline a drug-DESIGN tool, not just a screen. For each ranked lead it
predicts key ADMET endpoints from SMILES with ADMET-AI (Chemprop GNN, #1 on the TDC
leaderboard) and annotates the lead.

Deliberately ANNOTATION-ONLY: it adds fields, it does NOT change the ranking. Folding
ADMET into the composite score would re-open the validation-gate calibration; we
surface it as columns/flags first and can weight it later if we choose.

Auto-detects admet-ai and degrades gracefully (skips, labelled) if absent — same
discipline as the structure/pocket/docking nodes, so CI without the heavy model
still runs.
"""
from __future__ import annotations

from vta.admet.providers import ADMETAIProvider
from vta.state import VTAState

# raw ADMET-AI endpoint -> (friendly key, "higher better"?, risk threshold)
_ENDPOINTS = {
    "hERG": ("herg", False),                 # cardiotoxicity probability (lower better)
    "Bioavailability_Ma": ("oral", True),    # oral bioavailability 0-1 (higher better)
    "Solubility_AqSolDB": ("solubility", True),  # log mol/L (higher better)
}
_HERG_RISK = 0.5


def _admet_available() -> bool:
    return ADMETAIProvider().available()


def _row_value(preds, i, smiles, key):
    """Pull one endpoint from ADMET-AI output (DataFrame for a list, dict for a str)."""
    try:
        import pandas as pd
        if isinstance(preds, pd.DataFrame):
            return float(preds.iloc[i][key])
    except Exception:
        pass
    if isinstance(preds, dict):
        return float(preds.get(key))
    return None


def admet_node(state: VTAState) -> VTAState:
    leads = state.get("lead_candidates") or []
    if not leads:
        return state
    if not _admet_available():
        state["audit_trail"].append("[skip] ADMET: admet-ai not installed")
        return state

    provider = ADMETAIProvider()
    smiles = [l.get("smiles") or "" for l in leads]
    preds = provider.predict(smiles)

    flagged = 0
    for i, lead in enumerate(leads):
        admet = {k: preds[i].get(k) for k in ("herg", "oral", "solubility")}
        admet.update({
            "cyp_ddi_risk": None,
            "mitochondrial_toxicity": None,
            "dili": None,
            "ames_genotoxicity": None,
            "permeability": None,
            "microsomal_clearance": None,
            "plasma_protein_binding": None,
            "applicability_domain": preds[i].get("applicability_domain", "unknown"),
            "endpoint_status": "core_admet_ai_only; antiviral endpoints require provider",
        })
        herg = admet.get("herg")
        lead["admet"] = admet
        lead["admet_flag"] = "hERG risk" if (herg is not None and herg > _HERG_RISK) else ""
        if lead["admet_flag"]:
            flagged += 1

    top = leads[0]
    state["versions"]["admet"] = "ADMET-AI (Chemprop) + antiviral endpoint contract"
    state["audit_trail"].append(
        f"ADMET: annotated {len(leads)} leads ({flagged} hERG-flagged); "
        f"top lead {top.get('ligand')} hERG={top.get('admet', {}).get('herg')}, "
        f"oral={top.get('admet', {}).get('oral')}"
    )
    return state
