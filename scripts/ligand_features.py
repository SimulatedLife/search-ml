#!/usr/bin/env python
"""
Calculate the 25 ligand descriptors used by SEARCH-ML for every molecule in an
SDF file and write them to CSV (one row per molecule).

Usage
-----
    python ligand_features.py -i Ligand.sdf -o dataset.csv
    python ligand_features.py -i custom_molecules.sdf -o dataset.csv
    python ligand_features.py -i library.sdf -o dataset.csv --no-prepare   # legacy behaviour

Preparation (default on)
------------------------
The reference libraries were featurised from 3D, hydrogen-complete
structures. To keep user uploads consistent with them:

* molecules with 3D coordinates but no explicit hydrogens get hydrogens
  added at sensible positions;
* molecules that are 2D (all z = 0) or have no coordinates are protonated,
  embedded in 3D (ETKDGv3, fixed seed) and minimised (MMFF94, UFF fallback).

Molecules that cannot be sanitised or embedded are skipped and reported.
Exit code 6 if no molecule could be processed.
"""

import argparse
import csv
import os
import sys

import numpy as np

try:
    from tqdm import tqdm
except ImportError:  # tqdm is optional
    def tqdm(iterable=None, **_kwargs):
        return iterable

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, Descriptors3D

RDLogger.DisableLog("rdApp.warning")

RANDOM_SEED = 0xF00D

# Exact column order expected by the trained scaler / models
DESIRED_COLUMN_ORDER = [
    "PDB_File", "LigVol", "MolWt", "MolLogP", "TPSA", "NumHDonors",
    "NumHAcceptors", "NumRotatableBonds", "NumAromaticRings",
    "FractionCSP3", "NumHeteroatoms", "HeavyAtomCount", "NumAtoms",
    "MolMR", "BertzCT", "NumCarbons", "Asphericity", "Eccentricity",
    "InertialShapeFactor", "NPR1", "NPR2", "PMI1", "PMI2", "PMI3",
    "RadiusOfGyration", "SpherocityIndex",
]

DESCRIPTOR_FUNCS_2D = {
    "MolWt": Descriptors.MolWt,
    "MolLogP": Descriptors.MolLogP,
    "TPSA": Descriptors.TPSA,
    "NumHDonors": Descriptors.NumHDonors,
    "NumHAcceptors": Descriptors.NumHAcceptors,
    "NumRotatableBonds": Descriptors.NumRotatableBonds,
    "NumAromaticRings": Descriptors.NumAromaticRings,
    "FractionCSP3": Descriptors.FractionCSP3,
    "NumHeteroatoms": Descriptors.NumHeteroatoms,
    "BertzCT": Descriptors.BertzCT,
}

DESCRIPTOR_FUNCS_3D = {
    "LigVol": AllChem.ComputeMolVolume,
    "Asphericity": Descriptors3D.Asphericity,
    "Eccentricity": Descriptors3D.Eccentricity,
    "InertialShapeFactor": Descriptors3D.InertialShapeFactor,
    "NPR1": Descriptors3D.NPR1,
    "NPR2": Descriptors3D.NPR2,
    "PMI1": Descriptors3D.PMI1,
    "PMI2": Descriptors3D.PMI2,
    "PMI3": Descriptors3D.PMI3,
    "RadiusOfGyration": Descriptors3D.RadiusOfGyration,
    "SpherocityIndex": Descriptors3D.SpherocityIndex,
}


def count_carbons(mol):
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)


def molecule_name(mol, index):
    for prop in ("IDENTIFIER", "_Name"):
        if mol.HasProp(prop):
            name = " ".join(mol.GetProp(prop).split())
            if name:
                return name[:120]
    return f"Molecule_{index + 1}"


def has_3d_coordinates(mol):
    if mol.GetNumConformers() == 0:
        return False
    pos = mol.GetConformer().GetPositions()
    return len(pos) > 0 and float(np.abs(pos[:, 2]).max()) > 1e-3


def prepare_molecule(mol):
    """Return (molecule, action) with action in {'none', 'added_h', 'embedded'}."""
    has_h = any(a.GetAtomicNum() == 1 for a in mol.GetAtoms())
    is_3d = has_3d_coordinates(mol)
    if is_3d and has_h:
        return mol, "none"
    if is_3d:
        return Chem.AddHs(mol, addCoords=True), "added_h"
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = RANDOM_SEED
    if AllChem.EmbedMolecule(mol_h, params) == -1:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol_h, params) == -1:
            raise ValueError("3D embedding failed")
    try:
        if AllChem.MMFFHasAllMoleculeParams(mol_h):
            AllChem.MMFFOptimizeMolecule(mol_h, maxIters=2000)
        else:
            AllChem.UFFOptimizeMolecule(mol_h, maxIters=2000)
    except Exception:
        pass  # keep the embedded geometry
    return mol_h, "embedded"


def compute_features(mol, name):
    features = {"PDB_File": name}
    for key, func in DESCRIPTOR_FUNCS_2D.items():
        features[key] = func(mol)
    for key, func in DESCRIPTOR_FUNCS_3D.items():
        features[key] = func(mol)
    features["HeavyAtomCount"] = mol.GetNumHeavyAtoms()
    features["NumAtoms"] = mol.GetNumAtoms()
    features["MolMR"] = Crippen.MolMR(mol)
    features["NumCarbons"] = count_carbons(mol)
    return features


def calculate_features_sequentially(sdf_file_path, output_csv_file, batch_size=5000,
                                    prepare=True, max_molecules=0):
    """Featurise every molecule in ``sdf_file_path``. Returns (processed, failed)."""
    if not os.path.exists(sdf_file_path):
        print(f"Error: input SDF file not found: {sdf_file_path}")
        return 0, 0

    print(f"Opening SDF file: {sdf_file_path}")
    supplier = Chem.SDMolSupplier(sdf_file_path, sanitize=False, removeHs=False)
    batch, processed, failed = [], 0, 0
    actions = {"none": 0, "added_h": 0, "embedded": 0}
    failures = []

    with open(output_csv_file, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=DESIRED_COLUMN_ORDER, extrasaction="ignore")
        writer.writeheader()
        for i, mol in enumerate(tqdm(supplier, desc="Processing molecules")):
            if max_molecules and i >= max_molecules:
                print(f"Stopping after the first {max_molecules} molecules (limit).")
                break
            try:
                if mol is None:
                    raise ValueError("unreadable record")
                Chem.SanitizeMol(mol)
                name = molecule_name(mol, i)
                if prepare:
                    mol, action = prepare_molecule(mol)
                    actions[action] += 1
                batch.append(compute_features(mol, name))
                processed += 1
            except Exception as exc:
                failed += 1
                if len(failures) < 50:
                    failures.append(f"  #{i + 1}: {exc}")
                continue
            if len(batch) >= batch_size:
                writer.writerows(batch)
                batch = []
        if batch:
            writer.writerows(batch)

    print("\n--- Feature calculation complete ---")
    print(f"Molecules processed: {processed}")
    print(f"Molecules skipped:   {failed}")
    if prepare and processed:
        print(f"Preparation: {actions['none']} used as-is, {actions['added_h']} hydrogens added, "
              f"{actions['embedded']} embedded in 3D")
    if failures:
        print("Skipped molecules (first 50):")
        print("\n".join(failures))
    print(f"Features saved to '{output_csv_file}'")
    return processed, failed


def main():
    parser = argparse.ArgumentParser(
        description="Calculate SEARCH-ML ligand descriptors from an SDF file.")
    parser.add_argument("-i", "--input", required=True, help="Input SDF (one or many molecules)")
    parser.add_argument("-o", "--output", required=True, help="Output CSV")
    parser.add_argument("-b", "--batchsize", type=int, default=5000, help="Rows written per batch")
    parser.add_argument("--no-prepare", action="store_true",
                        help="Do not add hydrogens / embed 2D molecules (legacy behaviour)")
    parser.add_argument("--max-molecules", type=int, default=0, help="Process at most N molecules")
    args = parser.parse_args()

    processed, _failed = calculate_features_sequentially(
        args.input, args.output, args.batchsize, prepare=not args.no_prepare,
        max_molecules=args.max_molecules)
    if processed == 0:
        print("Error: no molecule could be processed. Check that the file is a valid SDF.")
        return 6
    return 0


if __name__ == "__main__":
    sys.exit(main())
