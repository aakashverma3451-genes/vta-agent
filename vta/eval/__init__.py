"""vta.eval — retrospective virtual-screening evaluation (SPEC #5).

Pure, dependency-free metrics for quantifying how well the pipeline separates known
actives from property-matched decoys: enrichment factor (EF), BEDROC, and ROC-AUC.
Kept free of heavy deps so they are fully hermetic-testable and reusable by both the
benchmark driver (`scripts/benchmark_enrichment.py`) and ad-hoc analysis.
"""
from vta.eval.metrics import bedroc, enrichment_factor, enrichment_report, roc_auc

__all__ = ["enrichment_factor", "roc_auc", "bedroc", "enrichment_report"]
