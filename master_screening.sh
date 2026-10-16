#!/usr/bin/env bash
# ==============================================================================
# SEARCH-ML : forward virtual screening (protein -> compound library)
# Author    : Dheeraj Kumar Chaurasia, SCFBio, IIT Delhi (dheeraj@scfbio-iitd.res.in)
# Version   : 2.0.0
#
# Usage:
#   ./master_screening.sh [options] <PDB_FILE> <LIGAND> <MODE>
#
#   PDB_FILE  protein-ligand complex in the current directory (e.g. 9kte.pdb)
#   LIGAND    reference ligand that defines the pocket:
#               9UM            first residue named 9UM
#               9UM:A:413      exact residue (RES:CHAIN:RESNUM)
#   MODE      DB | FDA | BIMP | <any key in config/databases.tsv> | SINGLE | CUSTOM
#               SINGLE needs single_molecule.{sdf,pdb,smi,mol} here (or -m FILE)
#               CUSTOM needs custom_molecules.sdf here (or -m FILE)
#
# Options:
#   -m FILE   molecule file for SINGLE / CUSTOM mode
#   -p FILE   property filter file ("Name: min,max" per line)
#   -n        no property filters (database modes use parameters/parameter.txt by default)
#   -o DIR    job directory to create (default: a new timestamped folder here)
#   -w N      parallel scoring processes (default: auto)
#   -k N      keep only the N best hits (default: all)
#   -h        show this help
#
# Examples:
#   ./master_screening.sh 9kte.pdb 9UM DB
#   ./master_screening.sh -n 9kte.pdb 9UM:A:413 FDA
#   ./master_screening.sh -m aspirin.smi 9kte.pdb 9UM SINGLE
#
# Outputs (in the job directory): results.csv, summary.json, job.log, error.log,
# PROT_LIG.pdb, leap.log and the COMPLETED / FAILED flag.
#
# Exit codes: 1 SEARCH_ML_HOME not set, 2 input missing, 3 protein < 25 residues,
#             4 SINGLE input missing, 5 CUSTOM input missing, 6 invalid mode,
#             see USAGE.md for the full list.
#
# Copyright (c) 2025 SCFBio, IIT Delhi. All rights reserved.
# ==============================================================================
set -o pipefail

usage() { sed -n '2,/^# =====.*$/{/^# ====/d;s/^# \{0,1\}//;p}' "$0" | sed -n '1,40p'; }

MOLS="" PARAMS="" NOFILTER=0 JOBDIR="" WORKERS="" TOPN=0
while getopts ":m:p:no:w:k:h" opt; do
  case "$opt" in
    m) MOLS="$OPTARG" ;;
    p) PARAMS="$OPTARG" ;;
    n) NOFILTER=1 ;;
    o) JOBDIR="$OPTARG" ;;
    w) WORKERS="$OPTARG" ;;
    k) TOPN="$OPTARG" ;;
    h) usage; exit 0 ;;
    :) echo "Error: option -$OPTARG needs a value." >&2; exit 6 ;;
    \?) echo "Error: unknown option -$OPTARG (use -h for help)." >&2; exit 6 ;;
  esac
done
shift $((OPTIND - 1))

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 [options] <PDB_FILE> <LIGAND> <MODE>   (use -h for help)" >&2
  exit 6
fi

if [ -z "${SEARCH_ML_HOME:-}" ] || [ ! -d "$SEARCH_ML_HOME" ]; then
  echo "Error: SEARCH_ML_HOME is not defined (or is not a directory)."
  echo "Set it once for the conda environment, then re-activate it:"
  echo "   conda env config vars set SEARCH_ML_HOME=/path/to/search-ml --name search-ml"
  echo "   conda deactivate && conda activate search-ml"
  exit 1
fi
echo "SEARCH_ML_HOME is set to: $SEARCH_ML_HOME"

PY="${SEARCHML_PYTHON:-python}"
command -v "$PY" >/dev/null 2>&1 || PY=python3

ARGS=(forward --pdb "$1" --ligand "$2" --mode "$3" --top "$TOPN")
[ -n "$MOLS" ] && ARGS+=(--molecules "$MOLS")
[ -n "$PARAMS" ] && ARGS+=(--params "$PARAMS")
[ "$NOFILTER" -eq 1 ] && ARGS+=(--no-filters)
[ -n "$JOBDIR" ] && ARGS+=(--job-dir "$JOBDIR")
[ -n "$WORKERS" ] && ARGS+=(--workers "$WORKERS")

date
"$PY" "$SEARCH_ML_HOME/scripts/pipeline.py" "${ARGS[@]}"
rc=$?
date
exit $rc
