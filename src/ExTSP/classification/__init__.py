"""
Classification utilities: ClinVar-based splits and bootstrap AUC for triplet-level scores.

Canonical column names: :mod:`ExTSP.classification.columns`.
"""

from ExTSP.classification import columns
from ExTSP.classification.auc_bootstrap import (
    BootstrapAUCResult,
    bootstrap_auc,
    bootstrap_auc_for_disease,
    build_positive_negative_triplets,
    run_one_bootstrap,
)


__all__ = [
    "BootstrapAUCResult",
    "bootstrap_auc",
    "bootstrap_auc_for_disease",
    "run_one_bootstrap",
]
