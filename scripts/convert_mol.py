#!/usr/bin/env python
"""
Convert a single small molecule (SDF / MOL / MOL2 / PDB / SMILES) into a
3D, hydrogen-complete ``Ligand.sdf`` ready for descriptor calculation.

Usage
-----
    python convert_mol.py <input> <output>
    python convert_mol.py query_ligand Ligand          # base name: tries .sdf .mol .mol2 .pdb .smi .smiles
    python convert_mol.py single_molecule.pdb Ligand.sdf
    python convert_mol.py --smiles "CC(=O)Oc1ccccc1C(=O)O" Ligand

The molecule is always re-embedded (RDKit ETKDGv3, fixed random seed for
reproducibility) and minimised with MMFF94 (UFF fallback), exactly like the
molecules in the pre-computed libraries.

Exit codes: 0 success, 2 input not found, 4 could not parse, 5 could not embed.
"""

import argparse
import os
import sys

EXTENSIONS = [".sdf", ".mol", ".mol2", ".pdb", ".smi", ".smiles", ".txt"]
RANDOM_SEED = 0xF00D


def find_input_file(name):
    """Return (path, ext) for an exact path or a base name + known extension."""
    if os.path.isfile(name):
        ext = os.path.splitext(name)[1].lower()
        return name, ext
    for ext in EXTENSIONS:
        candidate = name + ext
        if os.path.isfile(candidate):
            return candidate, ext
    return None, None


def read_molecule(path, ext):
    from rdkit import Chem

    if ext in (".smi", ".smiles", ".txt"):
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                parts = line.strip().split()
                if parts and not parts[0].startswith("#"):
                    mol = Chem.MolFromSmiles(parts[0])
                    if mol is not None and len(parts) > 1:
                        mol.SetProp("_Name", " ".join(parts[1:])[:80])
                    return mol
        return None
    if ext == ".pdb":
        return Chem.MolFromPDBFile(path, removeHs=False)
    if ext == ".mol2":
        return Chem.MolFromMol2File(path, removeHs=False)
    if ext == ".mol":
        return Chem.MolFromMolFile(path, removeHs=False)
    if ext == ".sdf":
        suppl = Chem.SDMolSupplier(path, removeHs=False)
        for mol in suppl:
            if mol is not None:
                return mol
        return None
    return None


def to_3d(mol):
    """Add hydrogens, embed in 3D and minimise. Returns the new molecule or None."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = RANDOM_SEED
    if AllChem.EmbedMolecule(mol, params) == -1:
        print("[!] Standard 3D embedding failed. Retrying with random coordinates...")
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) == -1:
            return None
    try:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
        else:
            raise ValueError("no MMFF parameters")
    except Exception:
        print("[!] MMFF94 not applicable. Falling back to UFF optimisation...")
        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
        except Exception as exc:  # keep the embedded geometry
            print(f"[!] UFF optimisation failed ({exc}); keeping the embedded geometry.")
    if name:
        mol.SetProp("_Name", name)
    return mol


def convert(input_name, output_file, smiles=None):
    try:
        from rdkit import Chem
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.warning")
    except ImportError:
        print("[!] RDKit is not installed in this Python environment (activate the search-ml env).")
        return 3

    if smiles:
        mol = Chem.MolFromSmiles(smiles.strip())
        source = "SMILES string"
    else:
        path, ext = find_input_file(input_name)
        if not path:
            print(f"[!] Error: no file found for '{input_name}' "
                  f"(looked for {', '.join(EXTENSIONS)})")
            return 2
        print(f"[*] Detected file: {path} (format: {ext})")
        try:
            mol = read_molecule(path, ext)
        except Exception as exc:
            print(f"[!] Parsing error: {exc}")
            mol = None
        source = path

    if mol is None:
        print(f"[!] Error: could not parse a molecule from {source}. "
              "Check that the structure / SMILES is valid.")
        return 4
    if mol.GetNumHeavyAtoms() < 2:
        print("[!] Error: the molecule has fewer than 2 heavy atoms.")
        return 4

    print(f"[*] Parsed {mol.GetNumHeavyAtoms()} heavy atoms. Generating 3D coordinates (ETKDGv3)...")
    mol3d = to_3d(mol)
    if mol3d is None:
        print("[!] Error: RDKit could not generate 3D coordinates for this molecule.")
        return 5

    if not output_file.lower().endswith(".sdf"):
        output_file += ".sdf"
    if not mol3d.HasProp("_Name") or not mol3d.GetProp("_Name").strip():
        mol3d.SetProp("_Name", "query_ligand")
    writer = Chem.SDWriter(output_file)
    writer.write(mol3d)
    writer.close()
    print(f"[+] Success! 3D structure saved as: {output_file} "
          f"(SMILES: {Chem.MolToSmiles(Chem.RemoveHs(mol3d))})")
    return 0


def main():
    parser = argparse.ArgumentParser(description="Convert one molecule to a 3D SDF for SEARCH-ML.")
    parser.add_argument("input_name", nargs="?", help="Input file, or base name without extension")
    parser.add_argument("output_name", help="Output SDF file (.sdf added if missing)")
    parser.add_argument("--smiles", help="Use this SMILES string instead of an input file")
    args = parser.parse_args()
    if not args.smiles and not args.input_name:
        parser.error("give an input file or --smiles")
    return convert(args.input_name, args.output_name, args.smiles)


if __name__ == "__main__":
    sys.exit(main())
