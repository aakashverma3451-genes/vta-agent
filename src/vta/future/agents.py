"""Abstract contracts for future VTA scientific platform agents."""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any


@dataclass(frozen=True)
class AgentRequest:
    query: str
    context: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class AgentResult:
    summary: str
    evidence: list[dict[str, Any]] = field(default_factory=list)
    confidence: float | None = None


class ScientificAgent(ABC):
    """Base contract for scientific reasoning agents."""

    name: str

    @abstractmethod
    def run(self, request: AgentRequest) -> AgentResult:
        """Execute the agent and return a structured result."""


class DrugKnowledgeGraphService(ABC):
    """Contract for querying and updating a drug discovery knowledge graph."""

    @abstractmethod
    def query(self, cypher: str, parameters: dict[str, Any] | None = None) -> list[dict]:
        """Run a graph query and return records."""

    @abstractmethod
    def upsert_evidence(self, records: list[dict[str, Any]]) -> int:
        """Insert or update evidence records and return the count written."""


class LiteratureAgent(ScientificAgent):
    """Finds and summarizes literature evidence for targets, drugs, and viruses."""


class TargetPrioritizationAgent(ScientificAgent):
    """Ranks candidate targets using biology, tractability, and evidence."""


class DrugRepurposingAgent(ScientificAgent):
    """Proposes repurposing candidates with mechanism and evidence traces."""


class ExperimentDesignAgent(ScientificAgent):
    """Suggests validation experiments and computational follow-up studies."""
