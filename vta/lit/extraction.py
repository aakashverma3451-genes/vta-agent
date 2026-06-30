"""Span-grounded literature relation extraction guards."""
from __future__ import annotations

NEGATION_MARKERS = (" not ", " no ", " failed to ", " did not ", " without ")


def normalize_entity(name: str, curie: str | None = None) -> dict:
    return {"name": name, "id": curie or f"UNRESOLVED:{name.lower().replace(' ', '_')}"}


def extract_relation(candidate: dict) -> dict | None:
    """Accept only grounded, non-negated candidate relations."""
    span = (candidate.get("evidence_span") or "").strip()
    if not span:
        return None
    text = f" {span.lower()} "
    if any(marker in text for marker in NEGATION_MARKERS):
        return None
    subj = normalize_entity(candidate.get("subject", ""), candidate.get("subject_id"))
    obj = normalize_entity(candidate.get("object", ""), candidate.get("object_id"))
    if subj["id"].startswith("UNRESOLVED") or obj["id"].startswith("UNRESOLVED"):
        return None
    return {
        "subject": subj,
        "type": candidate.get("type"),
        "object": obj,
        "source": candidate.get("source"),
        "confidence": float(candidate.get("confidence", 0.5)),
        "evidence_span": span,
        "evidence_tier": candidate.get("evidence_tier", "literature"),
        "preprint": bool(candidate.get("preprint", False)),
        "retraction_status": candidate.get("retraction_status", "unknown"),
    }
