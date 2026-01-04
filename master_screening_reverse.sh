#!/bin/bash

# ==============================================================================
# Script Name: master_reverse.sh (REVERSE SCREENING / SERIAL VERSION)
# Author: Dheeraj Kumar Chaurasia, SCFBio, IIT Delhi
# Email: dheeraj@scfbio-iitd.res.in
# Date: December 26, 2025
# Version: 1.1.0 (Target Identification / Homo Sapiens Only)
#
# Description:
#   This script performs reverse virtual screening to identify potential 
#   protein targets for a given small molecule.
#   
#   It screens the input molecule against the pre-calculated "HOMO" (Human)
#   protein dataset.
#   
#   USAGE GUIDELINES:
#   1. Ensure a small molecule file named 'query_ligand' (with ext .sdf, 
#      .pdb, or .smi) exists in the directory where you run this script.
#   2. Run this script directly. No arguments are needed.
#
#   Example: 
#      ./master_reverse.sh
#
# Copyright (c) 2025 Supercomputing Facility for Bioinformatics and Computational Biology (SCFBio), IIT Delhi.
# All rights reserved.
# ==============================================================================

date

# --- Assign variables ---
JOB_DIR=`pwd`
JOBID=`date +%s`

# ==============================================================================
# SECTION 1: JOB SETUP AND STATUS REPORTING
# ==============================================================================

# Exit immediately if a command fails.
set -e

# Check for Environment Variable
if [ -z "$SEARCH_ML_HOME" ]; then
    echo "Error: SEARCH_ML_HOME is not defined."
    echo "Please set the SEARCH-ML home path."
    exit 1
else
    echo "SEARCH_ML_HOME is set to: $SEARCH_ML_HOME"
fi

# Create the job directory
mkdir $JOBID

# CRITICAL: Navigate to the specific job directory.
cd "$JOB_DIR/$JOBID" || { echo "FATAL: Could not navigate to job directory $JOB_DIR/$JOBID"; exit 1; }

# --- Define a function to handle errors ---
handle_error() {
  local exit_code=$?
  echo "Job failed with exit code $exit_code at or near line $1." > error.log
  touch FAILED
  rm -f RUNNING
  exit $exit_code
}

# --- Set the trap ---
trap 'handle_error $LINENO' ERR

# --- Signal that the job has officially started ---
touch RUNNING
exec > >(tee -a job.log)
exec 2> >(tee -a job.log >&2)

echo "Job started at $(date)"
echo "Job Directory: $JOB_DIR/$JOBID"
echo "Mode: Reverse Screening (Target Identification)"
echo "Target Database: Homo Sapiens"

# ==============================================================================
# SECTION 2: INPUT PREPARATION
# ==============================================================================

# Define Paths
if [ -d "$SEARCH_ML_HOME" ]; then
    # Hardcoded to Homo Sapiens dataset
    export data_path="$SEARCH_ML_HOME/datasets/HOMO" 
    export scripts_path="$SEARCH_ML_HOME/scripts"
    export models="$SEARCH_ML_HOME/models"
else
    echo "Error: SEARCH-ML directory not found at $SEARCH_ML_HOME."
    exit 1
fi

# --- 1. Find Input Molecule in Parent Directory ---
echo "Searching for input molecule 'query_ligand'..."

INPUT_FILE=""
# Check parent directory (..) because we are now inside the JobID folder
if [ -f "../query_ligand.sdf" ]; then INPUT_FILE="../query_ligand.sdf"; 
elif [ -f "../query_ligand.pdb" ]; then INPUT_FILE="../query_ligand.pdb"; 
elif [ -f "../query_ligand.smi" ]; then INPUT_FILE="../query_ligand.smi"; 
fi

if [ -z "$INPUT_FILE" ]; then
    echo "Error: No valid input file found in parent directory ($JOB_DIR)."
    echo "Expected: query_ligand.sdf, query_ligand.pdb, or query_ligand.smi"
    exit 2
fi

echo "Found input: $INPUT_FILE"
cp "$INPUT_FILE" ./

# --- 2. Link Resources ---
echo "Linking scripts, models, and database..."
ln -s "$models"/* .
ln -s "$scripts_path"/* .
ln -s "$data_path"/* .

# ==============================================================================
# SECTION 3: CORE LOGIC
# ==============================================================================

echo "--- Step 1: Ligand Conversion ---"
# Convert the detected input file to a standardized 'Ligand.sdf'
python convert_mol.py query_ligand Ligand

echo "--- Step 2: Feature Calculation ---"
# Calculate physicochemical features for the ligand
python calculate_ligand_features.py -i Ligand.sdf -o dataset.csv

echo "--- Step 3: Reverse Screening Prediction ---"
# Run the reverse screening logic (Ligand vs Homo Proteins)
python reverse_screening.py > results.txt

# ==============================================================================
# SECTION 4: FINAL STATUS REPORTING
# ==============================================================================

echo "Job completed successfully at $(date)."
echo "Results saved to: $JOB_DIR/$JOBID/results.csv"

touch COMPLETED
rm -f RUNNING
date
exit 0
