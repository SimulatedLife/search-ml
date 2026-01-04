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
| `master_screening_reverse.sh` | Reverse screening / Target identification |

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

---

## Expected Directory Structure

SEARCH-ML must be installed using the following layout:

```text
search-ml/  
├── datasets/  
│   ├── drugbank/  
│   ├── fda/  
│   ├── bimp/  
│   └── HOMO/  
│
├── example/  
│   ├── 9kte.pdb  
│   ├── single_molecule.sdf  
│   ├── custom_molecules.sdf  
│   └── query_molecule.sdf 
│ 
├── scripts/  
│   ├── protein_features.py  
│   ├── pocket_features.py  
│   ├── merge_protein_features.py  
│   ├── ligand_features.py  
│   ├── screening.py  
│   ├── reverse_screening.py  
│   ├── convert_mol.py  
│   └── leap.cmd
│
├── models/  
│   ├── final_cat_model.joblib  
│   ├── final_lgbm_model.joblib  
│   ├── final_meta_model.joblib  
│   ├── final_rf_model.joblib  
│   ├── final_xgb_model.joblib  
│   └── scaler_final.joblib  
│
├── parameters/  
│   ├── final_ensemble_model_order.txt  
│   └── parameter.txt 
│ 
├── environment.yml
├── install.sh
├── README.md
├── USAGE.md
├── master_screening.sh
└── master_screening_reverse.sh
```

## Forward Screening Inputs

### Protein Input

- **Format:** PDB (`.pdb`)
- **Requirements:**
  - Must contain standard amino acid residues.
  - Protein chain length must be **≥ 25 amino acids**.
  - A co-crystallized ligand must be present in the structure (for pocket definition).

**Example:**
`1abc.pdb`

### Ligand Residue Code

- Passed as a command-line argument.
- Must match the specific 3-letter residue name of the ligand inside the PDB file.

**Example:**
`LIG`

---

## Ligand Input Files (Mode-Specific)

Depending on the chosen `MODE`, specific files may be required in the directory.

| Mode | Required File | Location | Notes |
| :--- | :--- | :--- | :--- |
| **DB** | *None* | Auto-linked | Uses internal DrugBank database |
| **FDA** | *None* | Auto-linked | Uses internal FDA database |
| **BIMP** | *None* | Auto-linked | Uses internal BIMP database |
| **SINGLE** | `single_molecule.sdf`<br>*(or .pdb / .smi)* | **Parent directory** | Must be renamed exactly |
| **CUSTOM** | `custom_molecules.sdf` | **Parent directory** | Must be a valid SDF file |
| **Reverse** | `query_ligand.sdf`<br>*(or .pdb / .smi)* | **Parent directory** | Used for Target Identification |

---

## Forward Virtual Screening (`master_screening.sh`)

### Syntax

```bash
./master_screening.sh <PDB_FILE> <LIGAND_CODE> <MODE>
```
### Argument Breakdown

The script accepts three mandatory positional arguments:

| Position | Argument | Description | Format / constraints | Example |
| :--- | :--- | :--- | :--- | :--- |
| **$1** | `PDB_FILE` | The input protein structure file. | Must be a valid `.pdb` file containing the receptor and a co-crystallized ligand. | `1abc.pdb` |
| **$2** | `LIGAND_CODE` | The Residue Name of the co-crystallized ligand. | 3-letter uppercase code matching the PDB entry. Used to define the active site. | `LIG`, `ATP`, `HEM` |
| **$3** | `MODE` | The screening library or input mode. | Must be one of: `DB`, `FDA`, `BIMP`, `SINGLE`, `CUSTOM`. | `DB` |

## Reverse Virtual Screening (`master_screening_reverse.sh`)

### Purpose

- Target identification.
- Ligand-based reverse screening.
- Screens against **Homo sapiens protein dataset only**.

---

### Syntax

```bash
./master_screening_reverse.sh
```

(No arguments required)

### Input Requirements

One of the following must exist in the parent directory:

* `query_ligand.sdf`
* `query_ligand.pdb`
* `query_ligand.smi`

### Workflow

1.  Ligand normalization → `Ligand.sdf`
2.  Ligand feature calculation
3.  Reverse screening prediction
4.  Target ranking

### Example

```bash
./master_screening_reverse.sh
```

## Output Files

| File | Description |
| :--- | :--- |
| `job.log` | Complete stdout and stderr. |
| `results.txt` | Screening predictions. |
| `dataset.csv` | Ligand feature dataset. |
| `args.txt` | Input arguments (forward screening only). |
| `error.log` | Error details. |

---

## Job State Flags

| File | Meaning |
| :--- | :--- |
| `RUNNING` | Job currently executing. |
| `COMPLETED` | Job finished successfully. |
| `FAILED` | Job terminated due to error. |

---

## Exit Codes

| Code | Meaning |
| :--- | :--- |
| **1** | `SEARCH_ML_HOME` not set. |
| **2** | Input file missing. |
| **3** | Protein length < 25 residues. |
| **4** | `SINGLE` mode input missing. |
| **5** | `CUSTOM` mode input missing. |
| **6** | Invalid screening mode. |

---

## Notes and Best Practices

- **Serial execution only:** Do not submit directly to batch schedulers.
- Ensure the **Ligand Residue Name** matches the PDB exactly.
- Do not rename internal scripts.
- Keep datasets and models read-only.
- **Always inspect `job.log` for diagnostics.**

---
