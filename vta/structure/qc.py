"""Structure quality-control contracts."""
from __future__ import annotations


def qc_structure(record: dict) -> dict:
    method = record.get("method")
    mean_plddt = record.get("mean_plddt")
    resolution = record.get("resolution")
    rfree = record.get("rfree")
    issues = []
    if mean_plddt is not None and float(mean_plddt) < 70:
        issues.append("low_plddt")
    if resolution is not None and float(resolution) > 3.5:
        issues.append("low_resolution")
    if rfree is not None and float(rfree) > 0.35:
        issues.append("high_rfree")
    if method in {"experimental_failed", "refused_too_long"} or not record.get("pdb_path"):
        issues.append("no_structure")
    return {
        "status": "pass" if not issues else "flagged",
        "issues": issues,
        "plddt": mean_plddt,
        "pae_status": "not_available",
        "molprobity_status": "not_available",
        "cryptic_pocket_status": "not_available",
        "ranking_active": False,
    }
