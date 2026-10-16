import pandas as pd
import numpy as np
import os
from Bio.PDB import PDBParser, is_aa
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import sys # Import sys to handle command-line arguments

# --- Configuration ---
BINDING_POCKET_RADIUS = 8.0 # In Angstroms.

# --- Self-contained three_to_one mapping to avoid import errors ---
d3to1 = {
    'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
    'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
}

# --- Helper Functions ---
def get_ligand_and_protein(structure, ligand_resname):
    """Separates a ligand and all protein residues from a Bio.PDB structure object."""
    ligand = None
    protein_residues = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname().strip() == ligand_resname:
                    ligand = residue
                elif is_aa(residue, standard=True):
                    protein_residues.append(residue)
    return ligand, protein_residues

def calculate_full_protein_features(protein_residues):
    """Calculates features for the entire protein from its list of residues."""
    if not protein_residues:
        return {}
        
    # Create a single sequence string from all residues using our local mapping
    try:
        protein_sequence = "".join([d3to1[res.get_resname()] for res in protein_residues])
    except KeyError as e:
        # Handle non-standard amino acids that are not in our dictionary
        print(f"Warning: Skipping non-standard residue {e} in full protein analysis.")
        protein_sequence = "".join([d3to1[res.get_resname()] for res in protein_residues if res.get_resname() in d3to1])

    if not protein_sequence:
        return {}

    # Use BioPython's ProteinAnalysis tool
    try:
        analysed_protein = ProteinAnalysis(protein_sequence)
    except Exception as e:
        print(f"Could not analyze protein sequence: {e}")
        return {}

    full_protein_features = {}
    full_protein_features['protein_total_residues'] = len(protein_sequence)
    full_protein_features['protein_molecular_weight'] = analysed_protein.molecular_weight()
    full_protein_features['protein_isoelectric_point'] = analysed_protein.isoelectric_point()
    full_protein_features['protein_aromaticity'] = analysed_protein.aromaticity()
    
    # Secondary structure fractions (estimated from amino acid content)
    helix, turn, sheet = analysed_protein.secondary_structure_fraction()
    full_protein_features['protein_helix_fraction'] = helix
    full_protein_features['protein_turn_fraction'] = turn
    full_protein_features['protein_sheet_fraction'] = sheet
    
    return full_protein_features

def calculate_pocket_features(protein_residues, ligand):
    """Calculates a comprehensive set of features for the protein binding pocket."""
    if not protein_residues or not ligand:
        return {}

    ligand_atoms = [atom for atom in ligand.get_atoms()]
    if not ligand_atoms: return {}
    ligand_center = sum(atom.coord for atom in ligand_atoms) / len(ligand_atoms)

    # Identify pocket residues
    pocket_residues = [res for res in protein_residues if any(np.linalg.norm(atom.coord - ligand_center) < BINDING_POCKET_RADIUS for atom in res.get_atoms())]

    if not pocket_residues:
        print(f"No pocket residues found within {BINDING_POCKET_RADIUS}Å of ligand.")
        return {}

    pocket_features = {}
    
    # Basic counts
    pocket_features['pocket_num_residues'] = len(pocket_residues)
    pocket_features['pocket_heavy_atoms'] = sum(1 for res in pocket_residues for atom in res.get_atoms() if atom.element != 'H')

    # Amino Acid Composition
    aa_composition = {'ALA':0,'CYS':0,'ASP':0,'GLU':0,'PHE':0,'GLY':0,'HIS':0,'ILE':0,'LYS':0,'LEU':0,'MET':0,'ASN':0,'PRO':0,'GLN':0,'ARG':0,'SER':0,'THR':0,'VAL':0,'TRP':0,'TYR':0}
    for residue in pocket_residues:
        res_name = residue.get_resname()
        if res_name in aa_composition: aa_composition[res_name] += 1
    
    for aa, count in aa_composition.items():
        pocket_features[f"pocket_{aa}_count"] = count

    # Physicochemical properties of the pocket
    polar_residues = ['SER', 'THR', 'CYS', 'ASN', 'GLN', 'TYR']
    nonpolar_residues = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'PRO', 'MET', 'PHE', 'TRP']
    aromatic_residues = ['PHE', 'TRP', 'TYR', 'HIS']
    charged_pos = ['LYS', 'ARG', 'HIS']
    charged_neg = ['ASP', 'GLU']
    
    pocket_features['pocket_polar_count'] = sum(1 for res in pocket_residues if res.get_resname() in polar_residues)
    pocket_features['pocket_nonpolar_count'] = sum(1 for res in pocket_residues if res.get_resname() in nonpolar_residues)
    pocket_features['pocket_aromatic_count'] = sum(1 for res in pocket_residues if res.get_resname() in aromatic_residues)
    pocket_features['pocket_charged_pos_count'] = sum(1 for res in pocket_residues if res.get_resname() in charged_pos)
    pocket_features['pocket_charged_neg_count'] = sum(1 for res in pocket_residues if res.get_resname() in charged_neg)
    pocket_features['pocket_net_charge'] = pocket_features['pocket_charged_pos_count'] - pocket_features['pocket_charged_neg_count']

    # Hydrophobicity
    kd_scale = {'ALA':1.8,'CYS':2.5,'ASP':-3.5,'GLU':-3.5,'PHE':2.8,'GLY':-0.4,'HIS':-3.2,'ILE':4.5,'LYS':-3.9,'LEU':3.8,'MET':1.9,'ASN':-3.5,'PRO':-1.6,'GLN':-3.5,'ARG':-4.5,'SER':-0.8,'THR':-0.7,'VAL':4.2,'TRP':-0.9,'TYR':-1.3}
    total_hydrophobicity = sum(kd_scale.get(res.get_resname(), 0.0) for res in pocket_residues)
    pocket_features['pocket_hydrophobicity_avg'] = total_hydrophobicity / len(pocket_residues) if pocket_residues else 0
    
    return pocket_features

def process_complex(pdb_file, ligand_resname):
    """Top-level function to process a complex and calculate all protein/pocket features."""
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("complex", pdb_file)
    except Exception as e:
        print(f"BioPython could not parse {pdb_file}: {e}")
        return None

    ligand, all_protein_residues = get_ligand_and_protein(structure, ligand_resname)
    
    if not ligand:
        print(f"Ligand '{ligand_resname}' not found in {pdb_file}")
        return None
    if not all_protein_residues:
        print(f"No protein residues found in {pdb_file}")
        return None

    # Calculate features for the pocket and the full protein
    pocket_feats = calculate_pocket_features(all_protein_residues, ligand)
    protein_feats = calculate_full_protein_features(all_protein_residues)
    if not pocket_feats:
        print(f"Error: no protein residues within {BINDING_POCKET_RADIUS} A of ligand '{ligand_resname}' "
              "- the selected group does not sit in a pocket.")
        return None
    if not protein_feats:
        print("Error: could not derive a protein sequence from the structure.")
        return None

    # Combine all features
    combined_features = {"PDB_File": os.path.basename(pdb_file)}
    combined_features.update(protein_feats)
    combined_features.update(pocket_feats)
    
    return combined_features

if __name__ == "__main__":

    # 2. Read the list of complexes
    if len(sys.argv) != 3:
        print("Usage: python protein_features.py <pdb_file_path> <ligand_resname>")
        print("Example: python protein_features.py PROT_LIG.pdb VWW")
        sys.exit(1) # Exit if the wrong number of arguments is provided

    pdb_file_arg = sys.argv[1]
    ligand_id_arg = sys.argv[2]

    # Create a list containing only the single complex from arguments
    complex_list = [(pdb_file_arg, ligand_id_arg)]

    # 3. Process all complexes
    all_protein_features = []
    if complex_list:
        print("--- Starting Expanded Protein & Pocket Feature Calculation ---")
        for pdb_file, ligand_id in complex_list:
            print(f"\nProcessing {pdb_file} for pocket around {ligand_id}...")
            features = process_complex(pdb_file, ligand_id)
            if features:
                all_protein_features.append(features)
            else:
                print(f"Failed to calculate features for {pdb_file}")

    # 4. Create and save the final DataFrame
    if all_protein_features:
        final_df = pd.DataFrame(all_protein_features)
        # Fill NaN values that might arise from failed calculations
        final_df.fillna(0, inplace=True)
        final_df.to_csv("protein_properties.csv", index=False)
        print("\n\n--- Success! ---")
        print("Expanded protein pocket feature calculation complete.")
        print(f"Saved {len(final_df)} entries to protein_properties.csv")
        print("\nShowing a sample of new columns:")
        sample_cols = ['PDB_File', 'protein_molecular_weight', 'protein_isoelectric_point', 'pocket_net_charge', 'pocket_aromatic_count']
        print(final_df[[c for c in sample_cols if c in final_df.columns]].head())
    else:
        print("\nNo protein features were generated.")
        sys.exit(10)

