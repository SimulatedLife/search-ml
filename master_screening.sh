#!/bin/bash

# ==============================================================================
# Script Name: master.sh (MASTER NODE / SERIAL VERSION)
# Author: Dheeraj Kumar Chaurasia, SCFBio, IIT Delhi
# Email: dheeraj@scfbio-iitd.res.in
# Date: December 26, 2025
# Version: 1.6.0 (Supports DRUGBANK/FDA/SINGLE/CUSTOM with Mol Processing)
#
# Description:
#   This script performs virtual screening on a provided PDB file. It supports
#   screening against standard databases (DrugBank, FDA) or user-provided 
#   molecule files (Single or Custom library).
#   
#   USAGE GUIDELINES:
#   1. This script is designed for SERIAL execution.
#   2. Run this directly on the master node or an interactive session.
#   3. Arguments:
#      $1 : PDB Filename (without extension)
#      $2 : Ligand Code (e.g., LIG)
#      $3 : Screening Mode ("DB", "FDA", "BIMP", "SINGLE", "CUSTOM")
#
#   File Requirements for Custom Modes:
#      - SINGLE: Requires 'single_molecule.{sdf|pdb|mol|smi}' in the directory.
#      - CUSTOM: Requires 'custom_molecules.sdf' in the directory.
#
#   Example: 
#      ./master.sh 1abc.pdb LIG DB
#      ./master.sh 1abc.pdb LIG DB        (Standard Screening)
#      ./master.sh 1abc.pdb LIG SINGLE    (Single Molecule Screening)
#
# Copyright (c) 2025 Supercomputing Facility for Bioinformatics and Computational Biology (SCFBio), IIT Delhi.
# All rights reserved.
# ==============================================================================

date

# --- Assign variables ---
JOB_DIR=`pwd`
JOBID=`date +%s`
PDB_FILE=$1
LIGAND_CODE=$2
MODE_ARG=$3      # Renamed for clarity (can be DB, FDA, SINGLE, CUSTOM)

# ==============================================================================
# SECTION 1: JOB SETUP AND STATUS REPORTING
# ==============================================================================

# Exit immediately if a command fails.
set -e

if [ -z "$SEARCH_ML_HOME" ]; then
    echo "Error: SEARCH_ML_HOME is not defined."
    echo "Please set the SEARCH-ML home path by executing the below commands:"
    echo "   conda env config vars set SEARCH_ML_HOME="$(pwd)" --name $ENV_NAME"
    echo "   conda deactivate $ENV_NAME"
    echo "   conda activate $ENV_NAME"
    echo "   conda env config vars set SEARCH_ML_HOME="$(pwd)" --name \$ENV_NAME"
    echo "   conda deactivate \$ENV_NAME"
    echo "   conda activate \$ENV_NAME"
    echo ""
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
echo "$1 , $2 , $3" >>args.txt
exec > >(tee -a job.log)
exec 2> >(tee -a job.log >&2)

echo "Job started at $(date)"
echo "Job Directory: $JOB_DIR/$JOBID"

# ==============================================================================
# SECTION 2: SCIENTIFIC LOGIC AND INPUT VALIDATION
# ==============================================================================

if [ -d "$SEARCH_ML_HOME" ]; then
    export data_path="$SEARCH_ML_HOME/datasets"
    export scripts_path="$SEARCH_ML_HOME/scripts"
    export models="$SEARCH_ML_HOME/models"
    export parameters="$SEARCH_ML_HOME/parameters"
else
    echo "Error: SEARCH-ML directory not found at $SEARCH_ML_HOME."
    exit 1
fi

# --- 1. Check PDB file existence ---
echo "Checking PDB file existence: $PDB_FILE"
if [ ! -f "../$PDB_FILE" ]; then
    echo "Error: PDB file '$PDB_FILE' not found in parent directory."
    exit 2
else
    echo "PDB file '$PDB_FILE' found."
    cp "../$PDB_FILE" ./
    
    # Protein Validation Logic
    aa_count=$(grep '^ATOM' "$PDB_FILE" | grep -E "ALA|ARG|ASH|ASN|ASP|CYM|CYS|CYX|GLH|GLN|GLU|GLH|GLY|HIS|HID|HIE|HIP|ILE|LEU|LYN|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL" | grep ' CA ' | cut -c 18-26 | sort -u | wc -l)
    if [ "$aa_count" -gt 24 ]; then
        echo "Protein length ($aa_count) is OK."
    else
        echo "Error: Protein size must be >= 25 amino acids. Found: $aa_count."
        exit 3
    fi
fi

# --- 2. Validate Mode and Database Selection ---
echo "Validating screening mode: $MODE_ARG"

SCREENING_TYPE=""    # Internal flag: "DATABASE" or "USER_FILE"
DATABASE_NAME=""     # Used if SCREENING_TYPE is DATABASE
INPUT_MOL_FILE=""    # Used if SCREENING_TYPE is USER_FILE

case "$MODE_ARG" in
  "DB") 
    SCREENING_TYPE="DATABASE"
    DATABASE_NAME="drugbank"
    echo "Mode Selected: Standard DrugBank Screening"
    ;;
  "FDA") 
    SCREENING_TYPE="DATABASE"
    DATABASE_NAME="fda"
    echo "Mode Selected: FDA Approved Drugs Screening"
    ;;
  "BIMP") 
    SCREENING_TYPE="DATABASE"
    DATABASE_NAME="bimp"
    echo "Mode Selected: BIMP Drugs Screening"
    ;;
  "SINGLE")
    SCREENING_TYPE="USER_FILE"
    echo "Mode Selected: Single Molecule Screening"
    
    # Autodetect format logic
    FOUND=0
    for ext in sdf pdb smi; do
        if [ -f "../single_molecule.$ext" ]; then
            echo "Found input file: single_molecule.$ext"
            cp "../single_molecule.$ext" ./
            INPUT_MOL_FILE="single_molecule.$ext"
            FOUND=1
            break
        fi
    done
    
    if [ $FOUND -eq 0 ]; then
        echo "Error: 'SINGLE' mode selected, but no valid file found."
        echo "Please ensure one of the following exists in the directory:"
        echo "  - single_molecule.sdf"
        echo "  - single_molecule.pdb"
        echo "  - single_molecule.smi"
        exit 4
    fi
    ;;
  "CUSTOM")
    SCREENING_TYPE="USER_FILE"
    echo "Mode Selected: Custom Library Screening"
    
    # Custom mode strict check
    if [ -f "../custom_molecules.sdf" ]; then
        echo "Found input file: custom_molecules.sdf"
        cp "../custom_molecules.sdf" ./
        INPUT_MOL_FILE="custom_molecules.sdf"
    else
        echo "Error: 'CUSTOM' mode selected, but 'custom_molecules.sdf' not found."
        echo "Custom mode currently only supports .sdf format."
        exit 5
    fi
    ;;
  *)
    echo "Error: Invalid argument '$MODE_ARG'."
    echo "Allowed options: DB, FDA, BIMP, SINGLE, CUSTOM"
    exit 6
    ;;
esac

# --- Define the core scientific workflow as a function ---
run_screening_logic() {
    echo "Preparing PDB and ligand files..."
    grep '^ATOM' "$PDB_FILE" | grep -E "ALA|ARG|ASH|ASN|ASP|CYM|CYS|CYX|GLH|GLN|GLU|GLH|GLY|HIS|HID|HIE|HIP|ILE|LEU|LYN|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL" > PROT.pdb
    grep -E "^ATOM|^HETATM" "$PDB_FILE" | grep "$LIGAND_CODE" >> PROT.pdb
    
    # Extract Ligand Residue Name
    LIG=$(echo "$LIGAND_CODE" | awk '{print $1}')
    echo "$LIG" >> lig

    echo "Loading Intel module and running tleap..."
    # module load intel/2018.3.222
    tleap -f leap.cmd
    # module purge
    echo "tleap finished."

    # --- PROTEIN FEATURE CALCULATIONS ---
    echo "Running Python scripts for protein feature calculation..."
    python protein_features.py PROT_LIG.pdb "$LIG"
    python pocket_features.py PROT_LIG.pdb "$LIG"
    echo "Feature calculation finished."

    echo "Merging Protein Parameters..."
    python merge_protein_features.py
    echo "Merging finished."

    # --- LIGAND/CUSTOM PROCESSING LOGIC ---
    # This block executes ONLY for SINGLE or CUSTOM modes
    if [ "$MODE_ARG" == "SINGLE" ]; then
        echo "Processing Single Molecule..."
        echo "Converting $INPUT_MOL_FILE to Ligand.sdf..."
        # Converts input (e.g., single_molecule.pdb) to Ligand.sdf
        python convert_mol.py "$INPUT_MOL_FILE" Ligand
        
        echo "Calculating features for Ligand.sdf..."
        python calculate_ligand_features.py -i Ligand.sdf -o dataset.csv

    elif [ "$MODE_ARG" == "CUSTOM" ]; then
        echo "Processing Custom Library..."
        echo "Calculating features for custom_molecules.sdf..."
        # Directly calculates features on the custom sdf file
        python calculate_ligand_features.py -i custom_molecules.sdf -o dataset.csv
    fi

    # --- FINAL SCREENING ---
    echo "Starting final virtual screening..."
    # The screening.py script uses the generated dataset.csv or linked DB files
    python screening.py > results.txt
    echo "Virtual screening finished."
}

# --- 3. Execute workflow (Serial Mode) ---
echo "Linking environment files..."

# Link models and scripts (Always required)
ln -s "$models"/* ./
ln -s "$scripts_path"/* ./

# Database & Parameters Linking Logic
if [ "$SCREENING_TYPE" == "DATABASE" ]; then
    echo "Database Mode: Linking database files and parameters."
    ln -s "$data_path/$DATABASE_NAME"/* ./
    ln -s "$parameters"/* ./
else
    echo "User File Mode: Skipping database and parameter linking."
    echo "Input molecule file: $INPUT_MOL_FILE"
fi

# Execute logic
run_screening_logic

# ==============================================================================
# SECTION 3: FINAL STATUS REPORTING
# ==============================================================================

echo "Job completed successfully at $(date)."

touch COMPLETED
rm -f RUNNING
date
exit 0
