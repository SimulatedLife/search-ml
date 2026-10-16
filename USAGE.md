# SEARCH-ML usage guide

**Version 2.0.0** · Author: Dheeraj Kumar Chaurasia, SCFBio, IIT Delhi (<dheeraj@scfbio-iitd.res.in>)

This guide covers the command-line tools. For installation see the
[README](README.md); for running web-portal jobs on a cluster see
[docs/HPC_SETUP.md](docs/HPC_SETUP.md).

- [1. Before you start](#1-before-you-start)
- [2. Forward screening](#2-forward-screening)
- [3. Reverse screening](#3-reverse-screening)
- [4. Input formats](#4-input-formats)
- [5. Property filters](#5-property-filters)
- [6. Outputs](#6-outputs)
- [7. Exit codes](#7-exit-codes)
- [8. Advanced use](#8-advanced-use)
- [9. Performance](#9-performance)

## 1. Before you start

```bash
conda activate search-ml
echo $SEARCH_ML_HOME                 # must print the repository path
python $SEARCH_ML_HOME/scripts/pipeline.py check
```

The tools are meant for a workstation, a login node or an interactive session.
To run many jobs on a cluster, either wrap them in your own batch script or use
the HPC runner in `hpc/`.

## 2. Forward screening

```text
master_screening.sh [options] <PDB_FILE> <LIGAND> <MODE>
```

Run it from the folder that contains the PDB file. A new job folder named after the
current Unix time (for example `1767001234/`) is created there.

### Arguments

| | Description | Example |
|---|---|---|
| `PDB_FILE` | Protein–ligand complex. Needs ATOM records for at least 25 standard amino-acid residues and the reference ligand (as HETATM or ATOM records). | `9kte.pdb` |
| `LIGAND` | Reference ligand that defines the pocket. `RES` selects the first residue with that name; `RES:CHAIN:RESNUM` selects one exact residue (use `_` for a blank chain, e.g. `LIG:_:1`). The 9-character PDB form `"VWW A 210"` is also accepted. | `9UM`, `9UM:A:413` |
| `MODE` | Ligand source (see below). | `FDA` |

| `MODE` | Screens against | Default filters |
|---|---|---|
| `DB` | DrugBank (10,041) | `parameters/parameter.txt` |
| `FDA` | FDA-approved drugs (2,066) | `parameters/parameter.txt` |
| `BIMP` | BIMP phytochemicals (105,909) | `parameters/parameter.txt` |
| any key in `config/databases.tsv` | that library, if installed (`chembl`, `coconut`, `molport`, `pubchem`, `zinc`, your own) | `parameters/parameter.txt` |
| `SINGLE` | one molecule (`-m FILE`, or `single_molecule.sdf/.mol/.mol2/.pdb/.smi` in the folder) | none |
| `CUSTOM` | your SDF library (`-m FILE.sdf`, or `custom_molecules.sdf` in the folder) | none |

### Options

| Option | Description |
|---|---|
| `-m FILE` | Molecule file for `SINGLE` / `CUSTOM`. |
| `-p FILE` | Property filter file (section 5). |
| `-n` | Do not apply property filters. |
| `-w N` | Number of parallel scoring processes. Default: 1 for small libraries, one per CPU above 200,000 molecules. |
| `-k N` | Keep only the N best-scoring molecules. |
| `-o DIR` | Job folder name instead of the time stamp. |
| `-h` | Help. |

### Examples

```bash
cd $SEARCH_ML_HOME/example
python ../scripts/prepare_complex.py --list 9kte.pdb        # see the ligands / ions
../master_screening.sh 9kte.pdb 9UM DB                      # DrugBank, default filters
../master_screening.sh -n 9kte.pdb 9UM:A:413 FDA            # FDA drugs, no filters
../master_screening.sh -m aspirin.smi 9kte.pdb 9UM SINGLE   # one molecule
../master_screening.sh -m library.sdf -k 1000 9kte.pdb 9UM CUSTOM
```

### What happens

1. **Prepare** (`prepare_complex.py`): keeps the protein residues (standard names and
   AMBER variants such as HID/HIE/CYX) and the chosen ligand residue, keeps the first
   alternate location of every atom, checks the 25-residue minimum, writes `PROT.pdb`.
2. **Standardise** (`tleap -f leap.cmd`, ff19SB/GAFF2): writes `PROT_LIG.pdb`.
3. **Protein and pocket descriptors** (`protein_features.py`, `pocket_features.py`):
   residue composition within 8 Å of the ligand centre, charge, hydrophobicity, pocket
   volume and SASA. Merged into `protein.csv`.
4. **Ligands:** library modes use the pre-computed `dataset.csv`; `SINGLE` converts the
   molecule to 3D (`convert_mol.py`) and `CUSTOM` prepares every SDF entry
   (`ligand_features.py`, section 4).
5. **Scoring** (`scoring.py`): filters, scales, predicts with the four base models and
   the meta-model, and writes `results.csv` sorted strongest first.

## 3. Reverse screening

```text
master_screening_reverse.sh [-q FILE | -s SMILES] [-k N] [-w N] [-o DIR] [LIBRARY]
```

| | Description |
|---|---|
| query | `query_ligand.sdf`, `.mol`, `.mol2`, `.pdb` or `.smi` in the current folder, or `-q FILE`, or `-s "SMILES"`. Exactly one molecule. |
| `LIBRARY` | Target library key from `config/libraries.tsv`. Default `HOMO` (*Homo sapiens*, 80,000 binding pockets). |

```bash
../master_screening_reverse.sh                                  # uses query_ligand.* here
../master_screening_reverse.sh -s "CC(=O)Oc1ccccc1C(=O)O" -k 200  # top 200 human targets
```

Steps: the molecule is converted to 3D (`convert_mol.py`), its 25 descriptors are
computed (`ligand_features.py`), and it is scored against every pocket in the library.
In `results.csv` the identifier `10GS_VWW_A_210.pdb` means PDB entry 10GS, pocket of
ligand VWW, chain A, residue 210.

About 3% of the HOMO pockets have no pocket-composition descriptors; as in version 1
they are scored with the models' built-in handling of missing values. Add
`--drop-incomplete` when calling `reverse_screening.py` directly to skip them.

## 4. Input formats

| Input | Accepted | Notes |
|---|---|---|
| Protein | `.pdb` | Only the first MODEL is used. Waters and other het groups are ignored. |
| Single molecule | `.sdf`, `.mol`, `.mol2`, `.pdb`, `.smi`/`.smiles` (first token = SMILES, rest = name) | Always re-embedded in 3D with RDKit ETKDGv3 (fixed seed, reproducible) and minimised with MMFF94 (UFF fallback). |
| Custom library | multi-molecule `.sdf` | Entries that already have 3D coordinates and hydrogens are used as they are; 3D entries without hydrogens get hydrogens added; 2D entries are protonated and embedded in 3D. Unreadable entries are skipped and listed in `job.log`. The molecule name comes from the `IDENTIFIER` property or the title line. |

The reference libraries were featurised from 3D, hydrogen-complete structures, which
is why user molecules are prepared the same way (`ligand_features.py --no-prepare`
restores the version 1 behaviour).

## 5. Property filters

A filter file has one property per line, `Name: min,max` (or `Name: value`):

```text
MolWt: 200,480
MolLogP: -0.4,5
TPSA: 0,131.6
NumHDonors: 0,5
```

Available properties: `MolWt`, `MolLogP`, `TPSA`, `NumHDonors`, `NumHAcceptors`,
`NumRotatableBonds`, `NumAromaticRings`, `NumHeteroatoms`, `HeavyAtomCount`,
`NumAtoms` (including hydrogens), `MolMR`, `NumCarbons`, `FractionCSP3`, `BertzCT`,
`LigVol`, `RadiusOfGyration`.

Filters are applied before prediction, so a strict filter also makes the run faster.
The default file `parameters/parameter.txt` combines the Lipinski, Ghose, Veber and Egan
ranges.

## 6. Outputs

| File | Description |
|---|---|
| `results.csv` | Ranked predictions (columns described in the README). With `-k N` only the best N rows. |
| `results.csv.gz` | Compressed copy, written when `results.csv` is larger than 50 MB. |
| `summary.json` | `rows_read`, `filtered_out`, `invalid_rows`, `rows_written`, filters, best hit, mean/std of the predictions, histogram of predicted pKi (bins of 0.1 from 0 to 15). |
| `job.json` | The job specification (inputs, mode, library, filters). |
| `state.json` | `RUNNING` / `COMPLETED` / `FAILED` / `CANCELLED`, timestamps, exit code. |
| `RUNNING`, `COMPLETED`, `FAILED` | Empty flag files (kept for compatibility with version 1 scripts). |
| `progress.json` | Current stage, percent, rows scored, ETA (updated every second while scoring). |
| `job.log` | Everything the steps printed. |
| `error.log` | On failure: what went wrong, a suggestion, and the last lines of the failing step. |
| `args.txt` | The command-line arguments (forward mode). |
| `PROT.pdb`, `PROT_LIG.pdb`, `leap.log`, `prepare.json` | Prepared complex. |
| `protein_properties.csv`, `pocket_properties.csv`, `protein.csv` | Protein/pocket descriptors. |
| `Ligand.sdf`, `dataset.csv` | Prepared user molecule(s) and their descriptors. |
| `parameter.txt` | The filters actually applied. |

## 7. Exit codes

| Code | Meaning |
|---:|---|
| 0 | Success |
| 1 | `SEARCH_ML_HOME` not set |
| 2 | Input file missing or unreadable |
| 3 | Protein shorter than 25 residues |
| 4 | `SINGLE` mode input missing |
| 5 | `CUSTOM` mode input missing |
| 6 | Invalid mode, library or option |
| 7 | Model files missing |
| 8 | Descriptor columns missing from a library |
| 9 | Out of memory |
| 10 / 11 / 12 | Protein descriptors / pocket descriptors / merging failed |
| 13 | Reference ligand not found |
| 14 | tleap failed |
| 15 | Molecule could not be parsed or embedded |
| 16 | No valid molecule in the input |
| 17 | Environment problem (RDKit, Biopython, tleap or ML libraries missing) |
| 18 | Library not installed on this machine |
| 20 | Internal error |
| 130 | Cancelled / stopped by the scheduler |

## 8. Advanced use

### The orchestrator

```bash
python $SEARCH_ML_HOME/scripts/pipeline.py forward --pdb X.pdb --ligand LIG --mode zinc --top 50000 --workers 16
python $SEARCH_ML_HOME/scripts/pipeline.py reverse --query mol.sdf --library homo
python $SEARCH_ML_HOME/scripts/pipeline.py run path/to/job_dir     # re-run a prepared job folder
python $SEARCH_ML_HOME/scripts/pipeline.py catalog                  # libraries and whether they are installed
python $SEARCH_ML_HOME/scripts/pipeline.py check                    # environment check
```

### The scoring engine alone

If you already have `protein.csv` (one row) and a ligand feature table:

```bash
python $SEARCH_ML_HOME/scripts/screening.py \
  --molecules features.csv --protein protein.csv --models $SEARCH_ML_HOME/models \
  --params my_filters.txt --output results.csv --summary summary.json --workers 8 --top 10000
```

`--molecules` accepts several files or a folder (its `parts/` sub-folder is used if
present; `.csv.gz` is read directly). Run `screening.py --help` for every option.

### Environment variables

| Variable | Effect |
|---|---|
| `SEARCH_ML_HOME` | Repository root (required by the shell scripts). |
| `SEARCHML_DATASETS` / `SEARCHML_MODELS` | Use libraries / models from another location. |
| `SEARCHML_DATABASES_TSV` / `SEARCHML_LIBRARIES_TSV` | Alternative catalog files. |
| `SEARCHML_CPUS` | CPUs the pipeline may use (otherwise taken from SLURM/PBS or the machine). |
| `SEARCHML_WORKERS` | Default number of scoring workers for `screening.py`. |
| `SEARCHML_GZIP_MB` | Size above which `results.csv.gz` is written (default 50). |
| `SEARCHML_PYTHON` | Python used by the shell scripts (default `python`). |

## 9. Performance

Scoring is the only step whose cost grows with the library; preparing the protein
takes a few seconds to a minute. Scoring times measured on the SCFBio cluster with the
library split over 32 parallel workers:

| Library | Molecules | Scoring time |
|---|---:|---:|
| DrugBank | 10,041 | 8 s |
| BIMP | 105,909 | 14 s |
| ZINC subset | 752,551 | 24 s |
| COCONUT | 695,114 | 24 s |
| ChEMBL | 2,068,782 | 36 s |
| MolPort | 5,733,129 | 1 min 13 s |
| PubChem | 95,310,581 | 15 min 12 s |

Each scoring worker needs roughly 2–3 GB of memory (the models are shared between
workers on Linux). Use `-k` for large libraries to keep `results.csv` manageable.
