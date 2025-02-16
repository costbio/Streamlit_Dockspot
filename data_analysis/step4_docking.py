import pandas as pd
from prody import *
import os
import tqdm
import glob
import subprocess
import seaborn as sns
import scipy.cluster.hierarchy as sch
import numpy as np
from openbabel import pybel
import re

# CONVERTING

def pdb_to_pdbqt(pdb_path, pdbqt_path, pH=7.4):
    """
    Convert a PDB file to a PDBQT file needed by docking programs of the AutoDock family.

    Parameters
    ----------
    pdb_path: str or pathlib.Path
        Path to input PDB file.
    pdbqt_path: str or pathlib.path
        Path to output PDBQT file.
    pH: float
        Protonation at given pH.
    """
    molecule = list(pybel.readfile("pdb", str(pdb_path)))[0]
    # add hydrogens at given pH
    molecule.OBMol.CorrectForPH(pH)
    molecule.addh()
    # add partial charges to each atom
    for atom in molecule.atoms:
        atom.OBAtom.GetPartialCharge()
    molecule.write("pdbqt", str(pdbqt_path), overwrite=True)
    return pdbqt_path


def calc_box(pdb_path, pocket_res_list):
    syst = parsePDB(pdb_path)
    pocket_res_list = pocket_res_list.split(' ')
    pocket_res_list = [el for el in pocket_res_list if len(el)>2]
    chains = np.unique([el[0] for el in pocket_res_list])
    selection_string = ""
    for chain in chains:
        resnums = [el[2:] for el in pocket_res_list if el.startswith(chain)]
        resnums_string = ' '.join(resnums)
        selection_string = selection_string+f"chain {chain} and resnum {resnums_string} or "

    selection_string = selection_string.rstrip(" or ")
    syst_box = syst.select(selection_string)
    box_coords = syst_box.getCoords()
    box_center = calcCenter(box_coords)
    box_min = np.min(box_coords,axis=0)
    box_max = np.max(box_coords,axis=0)
    return box_center, box_min, box_max


def run_smina(
    ligand_path, protein_path, out_path, pocket_center, pocket_size, smina_exe, num_poses=10, exhaustiveness=8
):
    """
    Perform docking with Smina.

    Parameters
    ----------
    ligand_path: str or pathlib.Path
        Path to ligand PDBQT file that should be docked.
    protein_path: str or pathlib.Path
        Path to protein PDBQT file that should be docked to.
    out_path: str or pathlib.Path
        Path to which docking poses should be saved, SDF or PDB format.
    pocket_center: iterable of float or int
        Coordinates defining the center of the binding site.
    pocket_size: iterable of float or int
        Lengths of edges defining the binding site.
    num_poses: int
        Maximum number of poses to generate.
    exhaustiveness: int
        Accuracy of docking calculations.

    Returns
    -------
    output_text: str
        The output of the Smina calculation.
    """
    try:
        result = subprocess.run(
            [
                smina_exe,
                "--ligand",
                str(ligand_path),
                "--receptor",
                str(protein_path),
                "--out",
                str(out_path),
                "--center_x",
                str(pocket_center[0]),
                "--center_y",
                str(pocket_center[1]),
                "--center_z",
                str(pocket_center[2]),
                "--size_x",
                str(pocket_size[0]),
                "--size_y",
                str(pocket_size[1]),
                "--size_z",
                str(pocket_size[2]),
                "--num_modes",
                str(num_poses),
                "--exhaustiveness",
                str(exhaustiveness),
                '--log', 'smina_debug.log'
            ],
            check=True,
            capture_output=True,
            text=True,  # needed to capture output text
        )
        return result.stdout, result.stderr
    except subprocess.CalledProcessError as e:
        print(f"Error running smina: {e}")
        print(f"Command: {e.cmd}")
        print(f"Return code: {e.returncode}")
        print(f"Output: {e.output}")
        print(f"Error output: {e.stderr}")
        raise


def parse_smina_log(text):
    # Extracting the relevant lines using regular expressions
    pattern = r'(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
    matches = re.findall(pattern, text)

    # Creating a DataFrame from the matches
    columns = ['mode', 'affinity (kcal/mol)', 'rmsd l.b.', 'rmsd u.b.']
    data = pd.DataFrame(matches, columns=columns)

    # Convert appropriate columns to numeric types
    data['mode'] = data['mode'].astype(int)
    data['affinity (kcal/mol)'] = data['affinity (kcal/mol)'].astype(float)
    data['rmsd l.b.'] = data['rmsd l.b.'].astype(float)
    data['rmsd u.b.'] = data['rmsd u.b.'].astype(float)

    # Return the DataFrame
    return data


def dock_ensemble( df_rep_pockets, ligand_folder, smina_exe, out_folder):
    
    list_outputs = list()

    for row in df_rep_pockets.iterrows():
        # Prepare receptor
        receptor_pdb = row[1]['File name']
        # Take only protein parts
        syst = parsePDB(receptor_pdb)
        protein = syst.select('protein')
        protein_pdb = os.path.join(out_folder,os.path.basename(receptor_pdb))
        writePDB(protein_pdb, protein)
        receptor_pdbqt = os.path.join(protein_pdb[:-4]+".pdbqt")
        pdb_to_pdbqt(protein_pdb, receptor_pdbqt)
        print(f"Prepared {receptor_pdbqt}")

        # Calculate box
        box_center, box_min, box_max = calc_box(protein_pdb, row[1]['residues'])
        box_size = [np.abs(box_max[0]-box_min[0]),np.abs(box_max[1]-box_min[1]),np.abs(box_max[2]-box_min[2])]
        
        # Run smina docking
        docking_output_folder = protein_pdb[:-4]+'_smina'
        if not os.path.exists(docking_output_folder):
            os.makedirs(docking_output_folder)

        ligands = glob.glob(ligand_folder+'/*.pdbqt')
        for lig_path in tqdm.tqdm(ligands):
            out_path = os.path.join(docking_output_folder,os.path.basename(lig_path)[:-6]+'_smina.sdf')
            output = run_smina(lig_path, receptor_pdbqt, out_path, box_center, box_size, smina_exe)
            df_output = parse_smina_log(output)
            df_output['library'] = ligand_folder
            df_output['ligand'] = os.path.basename(lig_path)[:-6]
            df_output['receptor'] = os.path.basename(receptor_pdb)
            df_output['receptor_path'] = receptor_pdbqt
            list_outputs.append(df_output)

    df_outputs = pd.concat(list_outputs,axis=0, ignore_index=True)
    return df_outputs










