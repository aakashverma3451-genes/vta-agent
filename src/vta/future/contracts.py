"""Shared request/response contracts for future VTA services."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass(frozen=True)
class EvidenceRecord:
    source: str
    claim: str
    entities: dict[str, str]
    score: float | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class PrioritizedItem:
    identifier: str
    label: str
    score: float
    rationale: str
    evidence: list[EvidenceRecord] = field(default_factory=list)
