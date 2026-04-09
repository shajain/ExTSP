"""
Run bootstrap AUC for one disease and setting from the command line::

    python -m ExTSP.classification --disease CM --setting 2 --n-bootstrap 50
"""
from __future__ import annotations

import argparse
import sys

from ExTSP.classification.auc_bootstrap import bootstrap_auc_for_disease
from ExTSP.config import Diseases


def _build_parser() -> argparse.ArgumentParser:
    disease_choices = sorted(Diseases.keys())
    p = argparse.ArgumentParser(
        description="Bootstrap AUC (IsoPath and exTSP) for variant–isoform–tissue triplets.",
    )
    p.add_argument(
        "--disease",
        required=True,
        choices=disease_choices,
        metavar="NAME",
        help=f"Disease sheet: {', '.join(disease_choices)}",
    )
    p.add_argument(
        "--setting",
        type=int,
        required=True,
        choices=[1, 2, 3, 4, 5],
        help="ClinVar / non-target configuration (1–5).",
    )
    p.add_argument(
        "--n-bootstrap",
        type=int,
        default=100,
        metavar="N",
        help="Number of bootstrap replicates (default: 100).",
    )
    p.add_argument(
        "--random-state",
        type=int,
        default=None,
        metavar="SEED",
        help="RNG seed for reproducibility.",
    )
    p.add_argument(
        "--balanced-negative",
        action="store_true",
        help="Balance benign negative variants to match pathogenic negative count when possible.",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    results, summary = bootstrap_auc_for_disease(
        args.disease,
        args.setting,
        balanced_negative=args.balanced_negative,
        n_bootstrap=args.n_bootstrap,
        random_state=args.random_state,
    )
    print(f"disease={args.disease!r} setting={args.setting} n_bootstrap={args.n_bootstrap}")
    print("summary:", summary)
    print(f"completed {len(results)} replicates")
    return 0


if __name__ == "__main__":
    sys.exit(main())
