"""VTA-Agent — Viral Target Assessment orchestrator (Module 2 + routing).

Wraps TaxonAgent (Module 1, the classifier) and drives an end-to-end pipeline:
classify → confidence-route → fold → pocket → dock → rank. Built contract-first
around `vta.state.VTAState`; the expensive biology (pockets, docking) is mocked
in Phase 1 and swapped for real tools in Phase 2.
"""
from __future__ import annotations

__version__ = "0.1.0"
