# 📖 Search-ML: Complete User Manual

**Version:** 1.3.1  
**Date:** December 14, 2025  
**Author:** Dheeraj Kumar Chaurasia (SCFBio, IIT Delhi)

---

## 1. Overview
**Search-ML** is a serial execution pipeline designed to automate the virtual screening of protein-ligand complexes against high-value databases (**DrugBank** and **FDA**). This tool integrates structure preparation (AmberTools), feature extraction (Python), and machine learning predictions into a single command.

This script is optimized for the **Master Node** or interactive workstations and does **not** require a job scheduler (like PBS/Slurm).

---

## 2. Syntax & Execution

### ⚠️ Critical Execution Rule
**Do NOT use `sh` to run this script.**
This script uses Bash-specific features (process substitution). Always use `bash` or execute directly.

* ✅ **Correct:** `bash master.sh ...`
* ✅ **Correct:** `./master.sh ...`
* ❌ **Incorrect:** `sh master.sh ...` (Will cause `syntax error near unexpected token >`)

### Command Structure
```bash
./master.sh <PDB_FILENAME> <LIGAND_CODE> <DATABASE>
```
## Step-by-Step Examples

This section provides concrete scenarios to help you understand how to run the script for different databases.

### Scenario A: Screening against DrugBank (DB)
**Goal:** You want to screen the protein-ligand complex `4dfr.pdb` against the DrugBank database. The ligand inside the PDB file is named `MTX`.

1.  **Check your files:**
    Ensure `4dfr.pdb` is in your current directory.
2.  **Run the command:**
    ```bash
    bash master.sh 4dfr.pdb MTX DB
    ```
3.  **Verify:**
    The script will create a directory (e.g., `1734177600`). Inside, check `results.csv` for the affinity scores.

### Scenario B: Screening against FDA Approved Drugs (FDA)
**Goal:** You have a target protein file named `target_prot.pdb` with a co-crystallized ligand named `LIG`. You want to find FDA-approved drugs that might bind to this pocket.

1.  **Check your files:**
    Ensure `target_prot.pdb` is in your current directory.
2.  **Run the command:**
    ```bash
    ./master.sh target_prot.pdb LIG FDA
    ```
    *(Note: We use `target_prot` with the `.pdb` extension)*
3.  **Verify:**
    Check the `job.log` inside the new folder to ensure `tleap` and feature calculation finished without errors.

### Scenario C: Handling Non-Standard File Names
**Goal:** Your file is named `my_experiment_v2.pdb` and the ligand is `DRG`.

1.  **Run the command:**
    ```bash
    bash master.sh my_experiment_v2.pdb DRG DB
    ```

## Output Explanation

When the job runs, it automatically creates a unique directory named after the **current Unix timestamp** (e.g., `1734177600`) to prevent overwriting previous results.

### Directory Structure & File Descriptions

Inside the job directory, you will find the following files:

| File | Description |
| :--- | :--- |
| **`results.csv`** | **The Main Output.** Contains the predicted binding affinity scores for your target against the selected database. |
| **`job.log`** | **Execution Log.** Captures all standard output (STDOUT) and errors (STDERR). Check this first if a job fails. |
| **`args.txt`** | **Input Record.** Logs the exact arguments (`PDB`, `LIG`, `DB`) used for this run. |
| **`PROT.pdb`** | **Cleaned Protein.** The protein structure extracted from your input file, stripped of the ligand and waters. |
| **`PROT_LIG.pdb`** | **Complex Structure.** The merged protein-ligand complex used for feature calculation. |
| **`leap.log`** | **AmberTools Log.** detailed output from `tleap`. Useful for debugging atom type or parameterization errors. |
| **`COMPLETED`** | **Success Flag.** An empty file created only if the script reaches the end successfully. |


## Workflow Logic

This is what happens under the hood when you execute the script:

1.  **Validation:**
    * Checks if the input PDB file exists.
    * Verifies that the protein chain is valid (> 24 amino acids).
2.  **Preparation:**
    * Splits the input PDB into `PROT.pdb` (Protein only) and `lig` (Ligand code).
    * Combines them into a clean `PROT_LIG.pdb`.
3.  **Parameterization:**
    * Runs `tleap` (AmberTools) to generate topology and check for missing atoms.
4.  **Featurization:**
    * Calculates physicochemical properties using RDKit/Python.
    * Calculates Pocket Volume and Solvent Accessible Surface Area (SASA).
5.  **Merging:**
    * Combines protein features and ligand features into a unified dataset.
6.  **Prediction:**
    * Runs the pre-trained Machine Learning models (`screening.py`) to generate affinity scores.
  
## Troubleshooting

### Error: `syntax error near unexpected token >`
* **Cause:** You ran the script with `sh scriptname.sh`.
* **Fix:** Run it with `bash scriptname.sh` or `./scriptname.sh`.

### Error: `Protein size must be >= 25 amino acids`
* **Cause:** The input PDB has fewer than 25 residues, or the script failed to detect standard amino acids (ALA, ARG, etc.).
* **Fix:** Check your PDB file format. Ensure it uses standard PDB naming conventions for residues.

### Error: `SEARCH_ML_HOME directory not found`
* **Cause:** The path in the script does not match your machine's file structure.
* **Fix:** Execute the following commands on the terminal:
* ```bash
  conda env config vars set SEARCH_ML="$(pwd)" --name search-ml
  conda activate search-ml
  ```
### Error: `HOME directory not found`
* **Cause:** The hardcoded path in the script does not match your machine's file structure.
* **Fix:** Open `master.sh` and edit the `export HOME="..."` line to point to your `search-ml` folder.

## Copyright & Contact

**Author:** Dheeraj Kumar Chaurasia  
**Affiliation:** Supercomputing Facility for Bioinformatics and Computational Biology (SCFBio), IIT Delhi  
**Email:** [dheeraj@scfbio-iitd.res.in](mailto:dheeraj@scfbio-iitd.res.in)

Copyright (c) 2025 SCFBio, IIT Delhi. All rights reserved.

