import os
import numpy as np
import re
from openbabel import pybel
import pandas as pd
from sklearn.preprocessing import MultiLabelBinarizer
from scipy.cluster import hierarchy
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objs as go

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

