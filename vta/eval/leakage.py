"""Leakage audit helpers for ML scorer promotion gates."""
from __future__ import annotations


def normalize_target_id(value: str | None) -> str:
    return (value or "").strip().upper()


def leakage_audit(training_targets: list[str], eval_targets: list[str]) -> dict:
    train = {normalize_target_id(t) for t in training_targets if normalize_target_id(t)}
    evals = {normalize_target_id(t) for t in eval_targets if normalize_target_id(t)}
    overlap = sorted(train & evals)
    return {
        "status": "fail_overlap" if overlap else "pass_no_overlap",
        "overlap": overlap,
        "n_training_targets": len(train),
        "n_eval_targets": len(evals),
        "promotable": not overlap,
    }


def gate_recalibration_criteria() -> dict:
    return {
        "required": [
            "property_matched_decoys",
            "scaffold_or_time_split",
            "leakage_audit_pass",
            "positive_control_gate_pass",
            "early_recognition_improves_without_le_regression",
        ],
        "primary_metrics": ["BEDROC(alpha=20)", "EF1%", "logAUC"],
        "ranking_policy": "DL terms remain annotation-only until all required gates pass.",
    }
