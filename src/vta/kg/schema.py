"""Initial VTA drug discovery knowledge graph schema definitions."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class GraphNode:
    label: str
    required: tuple[str, ...]
    optional: tuple[str, ...] = ()


@dataclass(frozen=True)
class GraphRelationship:
    type: str
    start: str
    end: str
    required: tuple[str, ...] = ()
    optional: tuple[str, ...] = ()


NODE_SCHEMAS: dict[str, GraphNode] = {
    "Virus": GraphNode("Virus", ("id", "name"), ("taxonomy", "genome_accession")),
    "Protein": GraphNode("Protein", ("id", "name"), ("sequence", "length_aa")),
    "Drug": GraphNode("Drug", ("id", "name"), ("smiles", "approval_status")),
    "Target": GraphNode("Target", ("id", "name"), ("target_class", "organism")),
    "Publication": GraphNode("Publication", ("id", "title"), ("doi", "year", "abstract")),
    "ClinicalTrial": GraphNode("ClinicalTrial", ("id", "title"), ("phase", "status")),
    "Pathway": GraphNode("Pathway", ("id", "name"), ("source",)),
}


RELATIONSHIP_SCHEMAS: dict[str, GraphRelationship] = {
    "TARGETS": GraphRelationship("TARGETS", "Drug", "Target", ("evidence",), ("score",)),
    "INHIBITS": GraphRelationship("INHIBITS", "Drug", "Protein", ("evidence",), ("ic50",)),
    "ASSOCIATED_WITH": GraphRelationship(
        "ASSOCIATED_WITH", "Target", "Virus", ("evidence",), ("score",)),
    "PARTICIPATES_IN": GraphRelationship("PARTICIPATES_IN", "Protein", "Pathway"),
    "SUPPORTED_BY": GraphRelationship(
        "SUPPORTED_BY", "Target", "Publication", ("claim",), ("confidence",)),
}
