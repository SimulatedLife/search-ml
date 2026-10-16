import pandas as pd
import numpy as np
import os
from Bio.PDB import PDBParser, is_aa
from scipy.spatial import KDTree
import sys # Import sys to handle command-line arguments

np.random.seed(42)
# --- Configuration ---
# The radius of the probe (e.g., a water molecule) to define the surface.
PROBE_RADIUS = 1.4 # Angstroms
# The number of random points to generate for the Monte Carlo estimation.
# Higher numbers increase accuracy but also computation time.
N_POINTS = 100000 # Increased for better surface area accuracy
# The maximum distance from the ligand's surface to consider as part of the pocket.
POCKET_RADIUS_CUTOFF = 6.0 # Angstroms
# The thickness of the shell used to estimate the surface area.
SASA_SHELL_THICKNESS = 0.1 # Angstroms

# A simple dictionary for van der Waals radii of common elements.
VDW_RADII = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
    'F': 1.47, 'P': 1.80, 'S': 1.80, 'CL': 1.75,
    'BR': 1.85, 'I': 1.98, 'HE': 1.40, 'LI': 1.82,
    'NA': 2.27, 'MG': 1.73, 'AL': 1.84, 'SI': 2.10,
    'K': 2.75, 'CA': 2.31, 'ZN': 1.39, 'CU': 1.40,
    'DEFAULT': 1.70 # Default for other elements
}

def calculate_pocket_properties(pdb_file, ligand_resname):
    """
    Calculates the active site volume and SASA using a probe-based Monte Carlo method.

    Args:
        pdb_file (str): Path to the complex PDB file.
        ligand_resname (str): The three-letter code of the ligand.

    Returns:
        dict: A dictionary containing the calculated volume and sasa, or None on failure.
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("complex", pdb_file)
    except Exception as e:
        print(f"BioPython could not parse {pdb_file}: {e}")
        return None

    # --- 1. Separate Ligand and Protein & Get Coordinates ---
    ligand, protein_atoms = None, []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname().strip() == ligand_resname:
                    ligand = residue
                elif is_aa(residue, standard=True):
                    protein_atoms.extend(atom for atom in residue.get_atoms() if atom.element.strip() != 'H')

    if not ligand or not protein_atoms:
        print(f"Could not find ligand '{ligand_resname}' or protein in {pdb_file}")
        return None

    ligand_atoms = list(ligand.get_atoms())
    if not ligand_atoms: return None
    
    protein_atom_coords = np.array([atom.coord for atom in protein_atoms])
    protein_atom_radii = np.array([VDW_RADII.get(atom.element.strip(), VDW_RADII['DEFAULT']) for atom in protein_atoms])
    
    # --- 2. Define Bounding Box ---
    ligand_coords = np.array([atom.coord for atom in ligand_atoms])
    min_coords = np.min(ligand_coords, axis=0) - POCKET_RADIUS_CUTOFF
    max_coords = np.max(ligand_coords, axis=0) + POCKET_RADIUS_CUTOFF
    box_dims = max_coords - min_coords
    box_volume = np.prod(box_dims)
    
    # --- 3. Generate Random Points in the Box (Monte Carlo) ---
    random_points = min_coords + np.random.rand(N_POINTS, 3) * box_dims

    # --- 4. Use KDTrees for Efficient Distance Searching ---
    protein_kdtree = KDTree(protein_atom_coords)
    ligand_kdtree = KDTree(ligand_coords)

    # --- 5. Filter Points to Define Pocket Volume and Surface ---
    
    # Calculate all distances first
    ligand_distances, _ = ligand_kdtree.query(random_points)
    protein_distances, closest_atom_indices = protein_kdtree.query(random_points)
    radii_of_closest_atoms = protein_atom_radii[closest_atom_indices]

    # Condition for a point to be inside the pocket (SAS definition)
    # Must be near the ligand AND not clash with the probe-expanded surface of the nearest atom.
    is_near_ligand_mask = ligand_distances < POCKET_RADIUS_CUTOFF
    is_outside_protein_mask = protein_distances > (radii_of_closest_atoms + PROBE_RADIUS)
    
    final_volume_mask = is_near_ligand_mask & is_outside_protein_mask
    final_volume_points_count = np.sum(final_volume_mask)
    
    # Condition for a point to be on the SASA surface shell
    # It must be inside the pocket, but only within a thin shell from the protein surface.
    is_in_shell_mask = protein_distances < (radii_of_closest_atoms + PROBE_RADIUS + SASA_SHELL_THICKNESS)
    final_sasa_mask = final_volume_mask & is_in_shell_mask
    final_sasa_points_count = np.sum(final_sasa_mask)

    # --- 6. Calculate Final Volume and SASA ---
    volume = (final_volume_points_count / N_POINTS) * box_volume
    sasa = (final_sasa_points_count / N_POINTS) * (box_volume / SASA_SHELL_THICKNESS)
    
    return {'pocket_volume': volume, 'pocket_sasa': sasa}

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python pocket_features.py <pdb_file_path> <ligand_resname>")
        print("Example: python pocket_features.py PROT_LIG.pdb VWW")
        sys.exit(1) # Exit if the wrong number of arguments is provided

    pdb_file_arg = sys.argv[1]
    ligand_id_arg = sys.argv[2]

    # Create a list containing only the single complex from arguments
    complex_list = [(pdb_file_arg, ligand_id_arg)]
    # --- END OF MODIFIED SECTION ---

    # 3. Process all complexes (now just the one from arguments)
    all_results = []
    if complex_list: # This condition will always be true now if arguments are provided
        print("--- Starting Active Site Volume & SASA Calculation ---")
        for pdb_file, ligand_id in complex_list: # Loop will run only once
            print(f"\nProcessing {pdb_file} for pocket around {ligand_id}...")
            properties = calculate_pocket_properties(pdb_file, ligand_id)
            if properties is not None:
                print(f"  > Calculated Volume: {properties['pocket_volume']:.2f} Å³")
                print(f"  > Calculated SASA:   {properties['pocket_sasa']:.2f} Å²")
                properties['PDB_File'] = pdb_file
                all_results.append(properties)
            else:
                print(f"  > Failed to calculate properties.")

    # 4. Create and save the final DataFrame
    if all_results:
        final_df = pd.DataFrame(all_results)
        final_df = final_df[['PDB_File', 'pocket_volume', 'pocket_sasa']] # Reorder columns
        # Name the output CSV based on the input PDB and ligand
        output_csv_name = "pocket_properties.csv"
        final_df.to_csv(output_csv_name, index=False)
        print("\n\n--- Success! ---")
        print("Pocket property calculation complete.")
        print(f"Saved {len(final_df)} entries to {output_csv_name}")
    else:
        print("\nNo properties were calculated.")
        sys.exit(11)
