# SEARCH-ML Virtual Screening Suite
**Usage Guide (`usage.md`)**

**Author:** Dheeraj Kumar Chaurasia  
**Affiliation:** Supercomputing Facility for Bioinformatics and Computational Biology (SCFBio), IIT Delhi  
**Email:** dheeraj@scfbio-iitd.res.in  

---

## Overview

SEARCH-ML is a **serial virtual screening framework** for structure-based drug discovery and target identification.  
It supports both **forward virtual screening** (protein → ligands) and **reverse virtual screening** (ligand → protein targets).

This repository provides two primary execution scripts:

| Script | Purpose |
|------|--------|
| `master_screening.sh` | Forward virtual screening (Protein–Ligand based) |
| `master_reverse_screening.sh` | Reverse screening / Target identification |

### Key Features

- Serial execution (master node / interactive session)
- Multiple screening modes:
  - DrugBank
  - FDA-approved drugs
  - BIMP (Bioactivity of Indian Medicinal Plants)
  - Single molecule
  - Custom ligand library
- Automatic job directory creation
- Robust error handling with job state flags
- Machine-learning–based prediction pipeline

---

## Intended Usage

- HPC master node
- Interactive compute session
- Local Linux/macOS workstation

**Not intended for direct batch submission unless externally wrapped.**


## Expected Directory Structure

SEARCH-ML must be installed using the following layout:

SEARCH_ML_HOME/  
├── datasets/  
│   ├── drugbank/  
│   ├── fda/  
│   ├── bimp/  
│   ├── HOMO/  
│   └── ...  
├── scripts/  
│   ├── protein_features.py  
│   ├── pocket_features.py  
│   ├── merge_protein_features.py  
│   ├── calculate_ligand_features.py  
│   ├── screening.py  
│   ├── reverse_screening.py  
│   ├── convert_mol.py  
│   └── ...  
├── models/  
│   ├── *.pkl  
│   ├── *.joblib  
│   └── ...  
├── parameters/  
│   ├── *.txt  
│   ├── *.csv  
│   └── ...  
├── master_screening.sh  
└── master_reverse_screening.sh  

All datasets, scripts, models, and parameters are dynamically linked at runtime.

---

## Job Execution Layout

Each execution creates a unique job directory using a UNIX timestamp:

working_directory/  
├── 1704029384/  
│   ├── job.log  
│   ├── results.txt  
│   ├── dataset.csv  
│   ├── RUNNING  
│   ├── COMPLETED / FAILED  
│   └── error.log  



## Environment Setup

### Mandatory Environment Variable

SEARCH-ML requires the following environment variable:

SEARCH_ML_HOME=/absolute/path/to/SEARCH_ML_HOME

Temporary example:

export SEARCH_ML_HOME=/home/user/SEARCH_ML_HOME

Recommended (Conda):

conda env config vars set SEARCH_ML_HOME=/home/user/SEARCH_ML_HOME --name <env_name>  
conda deactivate <env_name>  
conda activate <env_name>

---

## Execution Mode

- Serial execution only
- Intended for:
  - Master node execution
  - Interactive HPC jobs
  - Local Linux/macOS systems

---

## Software Requirements

### System Utilities

- bash
- grep
- awk
- cut
- sort
- tee

### Scientific Software

- Python ≥ 3.9
- RDKit
- NumPy
- Pandas
- scikit-learn
- AmberTools (`tleap`)

All Python utilities are executed from `$SEARCH_ML_HOME/scripts`.


## Forward Screening Inputs

### Protein Input

- Format: PDB
- Must contain:
  - Standard amino acid residues
  - At least **25 amino acids**
- Ligand must be present in the structure

Example:

1abc.pdb

---

### Ligand Residue Code

- Passed as command-line argument
- Must match residue name in PDB

Example:

LIG

---

## Ligand Input Files (Mode-Specific)

| Mode | Required File | Location |
|----|----|----|
| DB | None | Auto-linked |
| FDA | None | Auto-linked |
| BIMP | None | Auto-linked |
| SINGLE | single_molecule.sdf / .pdb / .smi | Parent directory |
| CUSTOM | custom_molecules.sdf | Parent directory |
| Reverse | query_ligand.sdf / .pdb / .smi | Parent directory |


## Forward Virtual Screening (`master_screening.sh`)

### Syntax

./master_screening.sh <PDB_FILE> <LIGAND_CODE> <MODE>

---

### Arguments

| Argument | Description |
|-------|------------|
| PDB_FILE | Protein PDB file (with extension) |
| LIGAND_CODE | Ligand residue name |
| MODE | Screening mode |

---

### Supported Screening Modes

| Mode | Description |
|----|------------|
| DB | DrugBank screening |
| FDA | FDA-approved drug screening |
| BIMP | BIMP compound screening |
| SINGLE | Single molecule screening |
| CUSTOM | Custom ligand library screening |


## Screening Mode Details

### DB / FDA / BIMP

- Uses precomputed ligand datasets
- Links ML models automatically
- Loads parameter files
- No ligand preprocessing required

---

### SINGLE Mode

- Auto-detects:
  - single_molecule.sdf
  - single_molecule.pdb
  - single_molecule.smi
- Converts input to Ligand.sdf
- Computes ligand features
- Screens against protein target

---

### CUSTOM Mode

- Requires custom_molecules.sdf
- Processes entire ligand library
- Computes descriptors for all molecules

---

### Forward Screening Examples

./master_screening.sh 1abc.pdb LIG DB  
./master_screening.sh 1abc.pdb LIG FDA  
./master_screening.sh 1abc.pdb LIG BIMP  
./master_screening.sh 1abc.pdb LIG SINGLE  
./master_screening.sh 1abc.pdb LIG CUSTOM  


## Reverse Virtual Screening (`master_reverse_screening.sh`)

### Purpose

- Target identification
- Ligand-based reverse screening
- Screens against **Homo sapiens protein dataset only**

---

### Syntax

./master_reverse_screening.sh

(No arguments required)

---

### Input Requirements

One of the following must exist in the parent directory:

query_ligand.sdf  
query_ligand.pdb  
query_ligand.smi  

---

### Workflow

1. Ligand normalization → Ligand.sdf  
2. Ligand feature calculation  
3. Reverse screening prediction  
4. Target ranking

---

### Example

./master_reverse_screening.sh



## Output Files

| File | Description |
|----|------------|
| job.log | Complete stdout and stderr |
| results.txt | Screening predictions |
| dataset.csv | Ligand feature dataset |
| args.txt | Input arguments (forward screening only) |
| error.log | Error details |

---

## Job State Flags

| File | Meaning |
|----|--------|
| RUNNING | Job currently executing |
| COMPLETED | Job finished successfully |
| FAILED | Job terminated due to error |

---

## Exit Codes

| Code | Meaning |
|----|--------|
| 1 | SEARCH_ML_HOME not set |
| 2 | Input file missing |
| 3 | Protein length < 25 residues |
| 4 | SINGLE mode input missing |
| 5 | CUSTOM mode input missing |
| 6 | Invalid screening mode |

---

## Notes and Best Practices

- Serial execution only
- Do not submit directly to batch schedulers
- Ensure ligand residue name matches PDB
- Do not rename internal scripts
- Keep datasets and models read-only
- Always inspect job.log for diagnostics

---

## Citation

If you use SEARCH-ML in your research, please cite:

Chaurasia DK et al.  
*Exploring chemical space for drug-like small molecules in the age of AI*  
Frontiers in Molecular Biosciences, 2025
