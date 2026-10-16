#!/usr/bin/env bash
# ==============================================================================
# SEARCH-ML : reverse virtual screening / target identification (molecule -> proteins)
# Author    : Dheeraj Kumar Chaurasia, SCFBio, IIT Delhi (dheeraj@scfbio-iitd.res.in)
# Version   : 2.0.0
#
# Usage:
#   ./master_screening_reverse.sh [options] [LIBRARY]
#
#   LIBRARY   target library key from config/libraries.tsv (default: HOMO = Homo sapiens)
#
# The query molecule is read from query_ligand.sdf / .pdb / .smi / .mol in the
# current directory, or given with -q FILE or -s SMILES.
#
# Options:
#   -q FILE    query molecule file
#   -s SMILES  query molecule as a SMILES string
#   -o DIR     job directory to create (default: a new timestamped folder here)
#   -w N       parallel scoring processes (default: auto)
#   -k N       keep only the N best targets (default: all)
#   -h         show this help
#
# Examples:
#   ./master_screening_reverse.sh
#   ./master_screening_reverse.sh -s "CC(=O)Oc1ccccc1C(=O)O" HOMO
#
# Output: results.csv (targets ranked by predicted affinity), summary.json, job.log.
# Copyright (c) 2025 SCFBio, IIT Delhi. All rights reserved.
# ==============================================================================
set -o pipefail

usage() { sed -n '2,/^# =====.*$/{/^# ====/d;s/^# \{0,1\}//;p}' "$0" | sed -n '1,30p'; }

QUERY="" SMILES="" JOBDIR="" WORKERS="" TOPN=0
while getopts ":q:s:o:w:k:h" opt; do
  case "$opt" in
    q) QUERY="$OPTARG" ;;
    s) SMILES="$OPTARG" ;;
    o) JOBDIR="$OPTARG" ;;
    w) WORKERS="$OPTARG" ;;
    k) TOPN="$OPTARG" ;;
    h) usage; exit 0 ;;
    :) echo "Error: option -$OPTARG needs a value." >&2; exit 6 ;;
    \?) echo "Error: unknown option -$OPTARG (use -h for help)." >&2; exit 6 ;;
  esac
done
shift $((OPTIND - 1))
LIBRARY="${1:-HOMO}"

if [ -z "${SEARCH_ML_HOME:-}" ] || [ ! -d "$SEARCH_ML_HOME" ]; then
  echo "Error: SEARCH_ML_HOME is not defined (or is not a directory)."
  echo "   conda env config vars set SEARCH_ML_HOME=/path/to/search-ml --name search-ml"
  echo "   conda deactivate && conda activate search-ml"
  exit 1
fi
echo "SEARCH_ML_HOME is set to: $SEARCH_ML_HOME"

PY="${SEARCHML_PYTHON:-python}"
command -v "$PY" >/dev/null 2>&1 || PY=python3

ARGS=(reverse --library "$LIBRARY" --top "$TOPN")
[ -n "$QUERY" ] && ARGS+=(--query "$QUERY")
[ -n "$SMILES" ] && ARGS+=(--smiles "$SMILES")
[ -n "$JOBDIR" ] && ARGS+=(--job-dir "$JOBDIR")
[ -n "$WORKERS" ] && ARGS+=(--workers "$WORKERS")

date
"$PY" "$SEARCH_ML_HOME/scripts/pipeline.py" "${ARGS[@]}"
rc=$?
date
exit $rc
