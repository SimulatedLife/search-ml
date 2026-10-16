#!/usr/bin/env python
"""
SEARCH-ML reverse screening (target identification): one molecule vs. a library
of pre-computed protein pockets.

Running it without arguments inside a job directory reads ``dataset.csv`` (the
query ligand features) and ``protein.csv`` (the target library) and writes
``results.csv``.

Examples
--------
    python reverse_screening.py
    python reverse_screening.py --ligand dataset.csv \\
        --proteins $SEARCH_ML_HOME/datasets/HOMO/protein.csv \\
        --models $SEARCH_ML_HOME/models --output results.csv --workers 8
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

import scoring  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(
        description="SEARCH-ML reverse screening (molecule vs. protein library)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--ligand", default="dataset.csv",
                   help="Query ligand features, one row (default: dataset.csv)")
    p.add_argument("--proteins", nargs="+", default=["protein.csv"],
                   help="Protein library feature CSV(s) or directory (default: protein.csv)")
    scoring.add_common_args(p)
    a = p.parse_args()

    def run():
        scoring.log("--- SEARCH-ML reverse virtual screening (1 ligand vs. many proteins) ---")
        scoring.run_screening(
            inputs=a.proteins,
            fixed_csv=a.ligand,
            models_dir=a.models,
            output=a.output,
            role="protein",
            params_file=None,
            order_file=a.order,
            preview=a.preview,
            preview_rows=a.preview_rows,
            summary=a.summary,
            progress=a.progress,
            workers=a.workers,
            threads=a.threads,
            top=a.top,
            chunksize=min(a.chunksize, 20_000),
            total_rows=a.total_rows,
            nan_policy="drop" if a.drop_incomplete else "auto",
        )
        scoring.log("--- Completed ---")

    return scoring.cli_main(run)


if __name__ == "__main__":
    sys.exit(main())
