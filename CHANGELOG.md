# Changelog

## 2.0.0 – 2026-10

### Fixed

- **SINGLE and CUSTOM modes failed every time**: `master_screening.sh` called
  `calculate_ligand_features.py`, which does not exist (the script is `ligand_features.py`).
- **Reverse screening failed**: `convert_mol.py` was called with the base name
  `query_ligand` and could not find the file; the final message pointed to
  `results.csv` while predictions were written elsewhere.
- **Meta-model input order** in SINGLE, CUSTOM and reverse modes: the order file was
  only linked in library modes, so the fallback order (`xgb, cat, …`) was used instead of
  the trained order (`cat, xgb, rf, lgbm`). Scores in those modes change slightly.
- Python steps printed errors but exited with status 0, so later steps failed with
  confusing messages; every step now exits non-zero with a clear reason.
- Ligand lines were selected with `grep LIGAND` anywhere in the line (could match atom
  names, other residues or both copies of a ligand); selection is now by PDB columns,
  with exact `RES:CHAIN:RESNUM` support. Only the first alternate location is kept.
- Job folders named after the current second collided when two runs started together.
- `install.sh`: environment name matched as a substring; `SEARCH_ML_HOME` unquoted;
  RDKit was missing from `environment.yml`.
- `merge_protein_features.py` silently wrote an empty `protein.csv` when the two
  feature files used different paths for `PDB_File`.
- `protein_features.py` crashed with a `KeyError` when no pocket residue was found;
  it now explains that the selected group is not in a pocket.

### Added

- `scripts/pipeline.py`: one orchestrator for the command line and the HPC runner
  (`forward`, `reverse`, `run`, `check`, `catalog`), with `state.json`, `progress.json`,
  `error.log` and documented exit codes.
- `scripts/scoring.py`: streaming, filter-before-predict, parallel workers, external
  merge sort, top-N limit, `.csv.gz` and split (`parts/`) libraries, summary statistics
  with an affinity histogram.
- `scripts/prepare_complex.py` (`--list` shows the ligands in a PDB file) and
  `scripts/dataset_tools.py` (count / split / validate / featurize libraries).
- User molecules are prepared like the reference libraries (hydrogens added, 2D
  entries embedded in 3D); 3D embedding uses a fixed seed (reproducible).
- `config/databases.tsv` and `config/libraries.tsv` library catalogs.
- `hpc/`: job manager for the web portal (SLURM, PBS, local), restricted SSH gateway,
  installer, self-test; `docs/HPC_SETUP.md`, `docs/ARCHITECTURE.md`.
- Unit and integration tests, GitHub Actions CI, `models/SHA256SUMS`, `CITATION.cff`.

### Changed

- `master_screening.sh` / `master_screening_reverse.sh` keep their positional
  arguments and exit codes 1–6, and gain options (`-m`, `-p`, `-n`, `-w`, `-k`, `-o`;
  reverse: `-q`, `-s`, library argument).
- Results are sorted by the unrounded predicted affinity (deterministic order for ties).
- `environment.yml` is now a short, cross-platform file; the machine-learning
  libraries stay pinned to the versions the models were saved with.

## 1.3.1 – 2025-12

- Initial public release of the serial pipeline.
