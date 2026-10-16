#!/usr/bin/env python
"""
SEARCH-ML forward virtual screening: one protein pocket vs. a molecule library.

Running it without arguments inside a job directory keeps the original
behaviour (reads ``dataset.csv`` + ``protein.csv`` + ``parameter.txt`` and the
models in ``./``, writes ``results.csv``).

Examples
--------
    python screening.py
    python screening.py --molecules $SEARCH_ML_HOME/datasets/fda/dataset.csv \\
        --protein protein.csv --models $SEARCH_ML_HOME/models \\
        --params parameter.txt --output results.csv --summary summary.json
    # large library split into parts, 16 worker processes, keep top 100k
    python screening.py --molecules /data/pubchem/parts --workers 16 --top 100000

Output columns: molecule ID, PBE (pKi/pKd/pIC50), PBE (kcal/mol), PBE (nM) and the
four base-model energies (XGB/CAT/RF/LGBM, kcal/mol), sorted strongest first.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

import scoring  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(
        description="SEARCH-ML forward screening (protein vs. molecule library)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--molecules", nargs="+", default=["dataset.csv"],
                   help="Ligand feature CSV(s) or a directory of CSV parts (default: dataset.csv)")
    p.add_argument("--protein", default="protein.csv", help="Protein + pocket features (default: protein.csv)")
    p.add_argument("--params", default="parameter.txt",
                   help="Property filter file 'Name: min,max' (default: parameter.txt; missing = no filters)")
    p.add_argument("--no-filters", action="store_true", help="Ignore any parameter file")
    scoring.add_common_args(p)
    a = p.parse_args()

    def run():
        scoring.log("--- SEARCH-ML forward virtual screening ---")
        scoring.run_screening(
            inputs=a.molecules,
            fixed_csv=a.protein,
            models_dir=a.models,
            output=a.output,
            role="ligand",
            params_file=None if a.no_filters else a.params,
            order_file=a.order,
            preview=a.preview,
            preview_rows=a.preview_rows,
            summary=a.summary,
            progress=a.progress,
            workers=a.workers,
            threads=a.threads,
            top=a.top,
            chunksize=a.chunksize,
            total_rows=a.total_rows,
            nan_policy="drop" if a.drop_incomplete else "auto",
        )
        scoring.log("--- Completed ---")

    return scoring.cli_main(run)


if __name__ == "__main__":
    sys.exit(main())
