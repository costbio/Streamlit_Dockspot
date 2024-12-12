import os
import numpy as np
import re
from openbabel import pybel
import pandas as pd
from sklearn.preprocessing import MultiLabelBinarizer
from scipy.cluster import hierarchy
import seaborn as sns
import matplotlib.pyplot as plt

# Define helper functions
def read_pdb(pdb_content):
    """
    Parse the PDB content to extract atomic coordinates.
    """
    coordinates = []
    for line in pdb_content.split('\n'):
        if line.startswith('ATOM'):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coordinates.append([x, y, z])
    return np.array(coordinates)

def center_of_mass(coordinates):
    """
    Calculate the center of mass for a set of coordinates.
    """
    total_mass = len(coordinates)
    com = np.sum(coordinates, axis=0) / total_mass
    return com

def create_grid_box(coordinates, padding=2.0):
    """
    Generate a grid box with specified padding around atomic coordinates.
    """
    min_coords = np.min(coordinates, axis=0)
    max_coords = np.max(coordinates, axis=0)

    grid_box_min = min_coords - padding
    grid_box_max = max_coords + padding

    return grid_box_min, grid_box_max

def keep_atoms(input_pdb_file, output_pdb_file, atom_numbers_to_keep):
    """
    Extract specific atoms from a PDB file based on atom numbers.
    """
    with open(input_pdb_file, 'r') as f:
        lines = f.readlines()

    with open(output_pdb_file, 'w') as f:
        for line in lines:
            if line.startswith('ATOM  '):
                atom_name = line[0:11].replace("ATOM", "").strip()
                if atom_name in atom_numbers_to_keep:
                    f.write(line)

def pdb_to_pdbqt(pdb_path, pdbqt_path, pH=7.4):
    """
    Convert a PDB file to PDBQT format for docking purposes.
    """
    molecule = list(pybel.readfile("pdb", str(pdb_path)))[0]
    for atom in molecule.atoms:
        atom.OBAtom.GetPartialCharge()
    molecule.write("pdbqt", str(pdbqt_path), overwrite=True)

def process_pdb_files(output_pdb_list, dock_locations):
    """
    Process a list of PDB files to calculate grid parameters and save them to text files.
    """
    for filename in output_pdb_list:
        try:
            with open(filename, 'r') as f:
                pdb_content = f.read()

            coordinates = read_pdb(pdb_content)
            com = center_of_mass(coordinates)

            respond_com = input("Use center of the pocket as center of mass? [yes/no]: ").lower()
            if respond_com == "yes":
                x_com, y_com, z_com = com
            else:
                x_com, y_com, z_com = map(float, input("Enter custom coordinates (x, y, z): ").split(','))

            respond_max = input("Use default grid box or custom? [yes/no]: ").lower()
            if respond_max == "yes":
                padding = float(input("Enter padding value (default 2.0): ") or 2.0)
                grid_box_min, grid_box_max = create_grid_box(coordinates, padding)

                max_min = input("Choose grid box [max/min]: ").lower()
                x_grid, y_grid, z_grid = (grid_box_max if max_min == "max" else grid_box_min)
            else:
                x_grid, y_grid, z_grid = map(float, input("Enter custom grid box coordinates (x, y, z): ").split(','))

            output_file_path = os.path.splitext(filename)[0] + '.txt'
            reseptor = dock_locations[output_pdb_list.index(filename)] + 'qt'

            with open(output_file_path, 'w') as output_file:
                output_file.writelines(f"receptor = {reseptor}\n")
                output_file.writelines(f"center_x = {x_com}\n")
                output_file.writelines(f"center_y = {y_com}\n")
                output_file.writelines(f"center_z = {z_com}\n")
                output_file.writelines(f"size_x = {x_grid}\n")
                output_file.writelines(f"size_y = {y_grid}\n")
                output_file.writelines(f"size_z = {z_grid}\n")
                output_file.writelines("energy_range = 3000\n")
                output_file.writelines("num_modes = 10\n")
                output_file.writelines("exhaustiveness = 8\n")

            print(f"Information saved to {output_file_path}")

        except FileNotFoundError:
            print(f"File {filename} not found.")

def load_data(csv_file_path):
    """
    Load the CSV file into a DataFrame.
    """
    return pd.read_csv(csv_file_path)

def preprocess_data(data):
    """
    Preprocess the data: split residues into lists and encode them.
    """
    residues = data['residues'].str.split()  # Split residues strings into lists
    mlb = MultiLabelBinarizer()
    residues_encoded = mlb.fit_transform(residues)
    return residues_encoded, mlb

def perform_clustering(residues_encoded, method='ward'):
    """
    Perform hierarchical clustering and return the linkage matrix.
    """
    linkage = hierarchy.linkage(residues_encoded, method=method)
    return linkage

def get_cluster_representatives(data, residues_encoded, linkage, depth):
    """
    Identify representative pockets for each cluster at the given depth.
    """
    cutree = hierarchy.cut_tree(linkage, n_clusters=depth)
    representatives = []

    def euclidean_distance(p, q):
        return np.sqrt(np.sum((p - q) ** 2))

    for cluster_id in range(depth):
        indices = [i for i, c in enumerate(cutree) if c[0] == cluster_id]
        if indices:
            cluster_points = residues_encoded[indices]
            centroid = np.mean(cluster_points, axis=0)
            distances = [euclidean_distance(centroid, point) for point in cluster_points]
            representative_index = np.argmin(distances)
            representative = cluster_points[representative_index]
            representative_pockets = data.iloc[indices[representative_index]]
            representatives.append((cluster_id, representative_pockets))

    return representatives

def save_representatives_to_file(representatives, file_path):
    """
    Save representative pockets to a text file.
    """
    with open(file_path, "w") as file:
        pass

    for rep in representatives:
        rep[1].to_csv(file_path, mode='a', index=True, header=False)

def cluster_and_save_representatives(csv_file_path, output_file_path, depth):
    """
    Full workflow: load data, preprocess, cluster, and save representatives.
    """
    data = load_data(csv_file_path)
    residues_encoded, _ = preprocess_data(data)
    linkage = perform_clustering(residues_encoded)
    representatives = get_cluster_representatives(data, residues_encoded, linkage, depth)
    save_representatives_to_file(representatives, output_file_path)
    return representatives

def main(processed_dir, mostsimilar_file_dir, csv_file_path, output_file_path, clustering_depth):
    """
    Main function to handle PDB file processing, docking preparation, and clustering analysis.
    """
    mostsimilar_file_dir = os.path.join(processed_dir, 'most_similar.txt')
    with open(mostsimilar_file_dir, "w") as file:
        content = file.read()

    blocks = content.strip().split('File name,')
    file_list = []
    residue_list = []
    data = {}

    for block in blocks[1:]:
        lines = block.strip().split('\n')
        file_name = lines[0]
        residues = lines[-1].split(", ")[1:]
        file_list.append(file_name)
        residue_list.append(residues)
        data[file_name] = residues

    pdb_list = ["{}.pdb".format(i) for i in file_list]
    pdb_files = f"{processed_dir}/4FOB_backbone_pdb_files"
    dock_locations = [f"{pdb_files}/{pdb}" for pdb in pdb_list]

    res_list = []
    for residues in residue_list:
        chosen = str(residues).replace("['", "").replace("']", "").replace("A_", "").split(" ")
        chosen = [element.replace('X_', '') for element in chosen]
        res_list.append(chosen)

    output_pdb_list = []
    for i in range(len(pdb_list)):
        output_pdb_file = f'output_{file_list[i]}.pdb'
        output_pdb_list.append(output_pdb_file)
        prot_loc = f"{pdb_files}/{pdb_list[i]}"
        keep_atoms(prot_loc, output_pdb_file, res_list[i])

    prot = "4FOB.pdb"
    output_pdb_list.append(prot)

    for a in dock_locations:
        pdb_to_pdbqt(a, a + "qt")

    process_pdb_files(output_pdb_list, dock_locations)

    # Perform clustering and save representatives
    cluster_and_save_representatives(csv_file_path, output_file_path, clustering_depth)

# Allow the script to be used as a module or standalone
if __name__ == "__main__":
    processed_dir = "path_to_processed_dir"
    mostsimilar_file_dir = f"{processed_dir}/most_similar.txt"
    csv_file_path = "path_to_csv_file"
    output_file_path = "path_to_output_file"
    clustering_depth = int(input("Enter clustering depth (default 3): ") or 3)
    main(processed_dir, mostsimilar_file_dir, csv_file_path, output_file_path, clustering_depth)
