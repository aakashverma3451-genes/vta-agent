"""Biolink-aligned VTA knowledge graph schema."""
from __future__ import annotations

NODES = {
    "Virus": "biolink:OrganismTaxon",
    "Protein": "biolink:Protein",
    "Gene": "biolink:Gene",
    "Variant": "biolink:SequenceVariant",
    "Drug": "biolink:Drug",
    "Target": "biolink:BiologicalEntity",
    "Publication": "biolink:Publication",
    "ClinicalTrial": "biolink:ClinicalTrial",
    "Pathway": "biolink:Pathway",
    "Disease": "biolink:Disease",
    "Assay": "biolink:ActivityAndBehavior",
    "Bioactivity": "biolink:Association",
    "Scaffold": "biolink:ChemicalEntity",
    "HostProtein": "biolink:Protein",
    "AdverseEvent": "biolink:DiseaseOrPhenotypicFeature",
    "Strain": "biolink:OrganismTaxon",
    "CellLine": "biolink:CellLine",
}

RELATIONSHIPS = {
    "TARGETS",
    "INHIBITS",
    "ASSOCIATED_WITH",
    "PARTICIPATES_IN",
    "SUPPORTED_BY",
    "CONFERS_RESISTANCE",
    "METABOLIZED_BY",
    "CAUSES_ADVERSE_EVENT",
    "EXPRESSED_IN",
    "INTERACTS_WITH",
    "HOST_FACTOR_FOR",
    "MEASURED_IN",
}

REQUIRED_EDGE_FIELDS = {"source", "confidence", "evidence_span"}

KG_SCHEMA = {
    "nodes": NODES,
    "relationships": sorted(RELATIONSHIPS),
    "edge_required_fields": sorted(REQUIRED_EDGE_FIELDS),
    "bootstrap_sources": ["Hetionet", "PrimeKG", "DRKG"],
    "identifier_policy": "CURIEs / identifiers.org",
}


def validate_edge(edge: dict) -> tuple[bool, list[str]]:
    missing = sorted(REQUIRED_EDGE_FIELDS - set(edge))
    if edge.get("type") not in RELATIONSHIPS:
        missing.append("valid_type")
    if not edge.get("evidence_span"):
        missing.append("non_empty_evidence_span")
    return not missing, missing
