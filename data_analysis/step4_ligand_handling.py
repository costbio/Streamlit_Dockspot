from openbabel import pybel
import tqdm
import os


def smiles_to_pdbqt(smiles, pdbqt_path, pH=7.4):
    """
    Convert a SMILES string to a PDBQT file needed by docking programs of the AutoDock family.

    Parameters
    ----------
    smiles: str
        SMILES string.
    pdbqt_path: str or pathlib.path
        Path to output PDBQT file.
    pH: float
        Protonation at given pH.
    """
    molecule = pybel.readstring("smi", smiles)
    # add hydrogens at given pH
    molecule.OBMol.CorrectForPH(pH)
    molecule.addh()
    # generate 3D coordinates
    molecule.make3D(forcefield="mmff94s", steps=10000)
    # add partial charges to each atom
    for atom in molecule.atoms:
        atom.OBAtom.GetPartialCharge()
    molecule.write("pdbqt", str(pdbqt_path), overwrite=True)
    return pdbqt_path



def prepare_ligands(df_ligands, out_folder):
    for row in tqdm.tqdm(df_ligands.iterrows()):
        smiles = row[1]['Ligand SMILES']
        bd_ln = row[1]['BindingDB MonomerID']
        if not os.path.exists(out_folder):
            os.makedirs(out_folder)
        pdbqt_path = os.path.join(out_folder,str(bd_ln)+'.pdbqt')
        smiles_to_pdbqt(smiles, pdbqt_path)


