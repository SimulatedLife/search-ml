#!/usr/bin/env python
"""
Merge protein_properties.csv (protein_features.py) and pocket_properties.csv
(pocket_features.py) into the single-row protein.csv used for screening.

Usage:  python merge_protein_features.py [protein_properties.csv pocket_properties.csv protein.csv]
Exit code 12 on failure.
"""

import sys

import pandas as pd


def merge_feature_csvs(protein_csv, pocket_csv, output_csv):
    """Merge the protein and pocket feature tables on 'PDB_File'. Returns True on success."""
    try:
        df_protein = pd.read_csv(protein_csv)
        print(f"Loaded '{protein_csv}' with shape: {df_protein.shape}")
    except FileNotFoundError:
        print(f"Error: the file '{protein_csv}' was not found (protein_features.py failed?).")
        return False

    try:
        df_pocket = pd.read_csv(pocket_csv)
        print(f"Loaded '{pocket_csv}' with shape: {df_pocket.shape}")
    except FileNotFoundError:
        print(f"Error: the file '{pocket_csv}' was not found (pocket_features.py failed?).")
        return False

    if len(df_protein) != 1:
        print("Error: protein feature CSV should contain exactly one row.")
        return False
    if len(df_pocket) != 1:
        print("Error: pocket feature CSV should contain exactly one row.")
        return False

    print("Merging protein features on 'PDB_File'...")
    df_merged = pd.merge(df_protein, df_pocket, on="PDB_File", how="inner")
    if df_merged.empty:
        # both files describe the same single complex; fall back to a column-wise join
        df_merged = pd.concat([df_protein.reset_index(drop=True),
                               df_pocket.drop(columns=["PDB_File"]).reset_index(drop=True)], axis=1)
    print(f"Shape after merge: {df_merged.shape}")

    df_merged.to_csv(output_csv, index=False)
    print("\n--- Success! ---")
    print(f"Merged all features into '{output_csv}'")
    return True


if __name__ == "__main__":
    args = sys.argv[1:]
    if args and len(args) != 3:
        print(__doc__)
        sys.exit(1)
    protein_features_file, pocket_features_file, final_output_file = (
        args if args else ["protein_properties.csv", "pocket_properties.csv", "protein.csv"])
    ok = merge_feature_csvs(protein_features_file, pocket_features_file, final_output_file)
    sys.exit(0 if ok else 12)
