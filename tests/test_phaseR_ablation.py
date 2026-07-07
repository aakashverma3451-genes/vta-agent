"""Tests for the Phase R / R6 component ablation. Hermetic — asserts routing-type-dependence.

Drives the real node logic via `_disposition` (reads committed gate artifacts, writes nothing),
so it proves the shipped behaviour: the router fixes out-of-domain targets, the verification
gate fixes in-domain-but-loses-to-2D targets, and neither alone is sufficient.
"""
from __future__ import annotations

from scripts.phaseR_ablation import TARGETS, _disposition

_BY_NAME = {t["name"]: t for t in TARGETS}


def _d(level, name):
    return _disposition(level, _BY_NAME[name])


def test_raw_pipeline_always_claims_enrichment():
    for name in _BY_NAME:
        assert _d("L0_raw", name) == "dock_and_claim_enrichment"


def test_router_fixes_out_of_domain_targets_at_L2():
    # the router turns unwarranted docking claims into honest annotate/defer/refuse
    assert _d("L2_triage", "NS5B") == "annotate_only"   # un_benchmarkable nucleotide
    assert _d("L2_triage", "GPX") == "refuse"           # binding-site pLDDT 45 < 50
    assert _d("L2_triage", "GPY") == "defer"            # binding-site pLDDT 60 borderline
    # ...but an in-domain target still claims enrichment at L2 (no gate yet)
    assert _d("L2_triage", "MPRO") == "dock_and_claim_enrichment"


def test_verification_fixes_in_domain_losers_at_L3():
    # the hard gate downgrades docking claims that don't beat the 2D baseline / lack a pose
    assert _d("L3_verification", "MPRO") == "dock_but_downgraded"
    assert _d("L3_verification", "PB1") == "dock_but_downgraded"


def test_every_target_correct_only_at_full_agent():
    for t in TARGETS:
        assert _d("L4_full_agent", t["name"]) == t["truth"], t["name"]


def test_no_single_layer_is_sufficient():
    # router alone leaves Mpro/PB1 over-claiming; verification alone can't route NS5B/GPX/GPY.
    l2 = {n: _d("L2_triage", n) for n in _BY_NAME}
    l2_correct = sum(l2[t["name"]] == t["truth"] for t in TARGETS)
    l4_correct = sum(_d("L4_full_agent", t["name"]) == t["truth"] for t in TARGETS)
    assert 0 < l2_correct < l4_correct == len(TARGETS)   # strict gain, neither layer alone enough
