import streamlit as st
import os
import json
import pandas as pd
import stmol
from data_analysis.step4_docking import *
from data_analysis.step4_ligand_handling import smiles_to_pdbqt
import time
import py3Dmol
import numpy as np
import subprocess

if 'new_job_id' not in st.session_state:
    st.session_state.new_job_id = None

# Base folder for jobs
BASE_JOB_FOLDER = 'streamlit_jobs'

STATUS_FOLDER = 'job_statuses'

def get_status_file(job_id):
    return os.path.join(STATUS_FOLDER, f'{job_id}.json')

def load_status(job_id):
    status_file = get_status_file(job_id)
    if os.path.exists(status_file):
        with open(status_file, 'r') as f:
            return json.load(f)
    return {}

def update_status(newjob_id, status, step=None):
    statuses = load_status(newjob_id)
    statuses['status'] = status
    statuses['step'] = step
    
    with open(get_status_file(newjob_id), 'w') as f:
        json.dump(statuses, f)

def get_jobid(jobid):
    file = os.path.join(STATUS_FOLDER, f'{jobid}.json')
    
    with open(file, 'r') as file:
        data = json.load(file)
    job_id = data["first_job_id"]
    return job_id




# Function to clean the pdbqt_files folder
def clean_pdbqt_folder(folder_path):
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if os.path.isfile(file_path):
            os.remove(file_path)

# Streamlit UI
#st.set_page_config(page_title="Docking App", layout="wide")
st.title("Docking App")


st.header("Job Configuration")
job_id_input = st.text_input("Enter your Job ID", value=st.session_state.get("new_job_id", ""))
smiles_input = st.text_area("Enter SMILES strings for ligands (one per line):", 
                                help="Enter SMILES strings for the ligands, separated by lines.")
run_button = st.button("Run Ensemble Docking")

# Main content
if job_id_input:
    if st.session_state.new_job_id is None:
        st.session_state.new_job_id = job_id_input
    job_folder = os.path.join(BASE_JOB_FOLDER, job_id_input)
    old_job_folder = os.path.join(BASE_JOB_FOLDER, get_jobid(st.session_state.new_job_id))
    if os.path.exists(job_folder):
        st.success(f"Job found: {job_id_input}")

        # Display job status
        job_status = load_status(job_id_input).get('status', 'Not Started')
        st.write(f"Job Status: {job_status}")
        oldjobid = get_jobid(st.session_state.new_job_id)
        st.write(oldjobid)

        if run_button:
            if not smiles_input:
                st.warning("Please enter SMILES strings to proceed.")
            else:
                ligands_smiles = [smiles.strip() for smiles in smiles_input.split('\n')]
                df = pd.DataFrame(ligands_smiles, columns=["Ligands"])
                st.dataframe(df, hide_index=True, use_container_width=True)
                st.info("Converting ligands...")

                # Creating necessary folders
                processed_dir = os.path.join(job_folder, "processed")
                old_processed_dir = os.path.join(old_job_folder, "processed")
                os.makedirs(processed_dir, exist_ok=True)

                pdbqt_folder = os.path.join(processed_dir, "pdbqt_files")
                os.makedirs(pdbqt_folder, exist_ok=True)

                # Clean previous pdbqt files
                clean_pdbqt_folder(pdbqt_folder)
                ligand_folder = os.path.join(pdbqt_folder, "ligand_pdbqt_files")
                os.makedirs(ligand_folder, exist_ok=True)

                prot_pdbqt_folder = os.path.join(pdbqt_folder,"protein_pdbqt_files")
                os.makedirs(prot_pdbqt_folder, exist_ok=True)

                # Loop through SMILES and process each ligand
                progress_bar = st.progress(0)
                for i, smiles in enumerate(ligands_smiles):
                    pdbqt_path = os.path.join(ligand_folder, f"ligand_{smiles}.pdbqt")
                    result_message = smiles_to_pdbqt(smiles, pdbqt_path)
                    progress_bar.progress((i + 1) / len(ligands_smiles))
                    time.sleep(1)  # simulate processing time

                st.success("All ligands are processed successfully!")

                # Step 3: Docking Setup
                st.header("Ensemble Docking")
                # Load the chosen representatives dataframe
                csv_path = os.path.join(processed_dir, "chosen_representatives.csv")
                df = pd.read_csv(csv_path)
                new_df = df[["Frame","pocket_index","probability","residues"]]
                st.write("Chosen receptors")
                st.dataframe(new_df, hide_index=True, use_container_width=True)
                
                out_folder = os.path.join(processed_dir, "docking_results")
                os.makedirs(out_folder, exist_ok=True)

                # Iterate through receptors and ligands
                list_outputs = []
                total_simulations = len(ligands_smiles) * len(df)
                completed_simulations = 0
                docking_progress_bar = st.progress(0)

                for index, row in df.iterrows():
                    receptor_pdb = os.path.join(old_processed_dir, "pdb_files", f"{row['File name']}.pdb")
                    from prody import *
                    syst = parsePDB(receptor_pdb)
                    protein = syst.select('protein')

                    if protein is None:
                        st.warning(f"No protein parts found in {receptor_pdb}")
                        continue
                    
                    protein_pdb = os.path.join(out_folder, os.path.basename(receptor_pdb))
                    os.makedirs(os.path.dirname(protein_pdb), exist_ok=True)
                    writePDB(protein_pdb, protein)

                    receptor_pdbqt = protein_pdb[:-4] + ".pdbqt"
                    receptor_pdbqt = receptor_pdbqt.replace("docking_results", "pdbqt_files/protein_pdbqt_files")
                    pdb_to_pdbqt(protein_pdb, receptor_pdbqt)

                    filename = os.path.basename(receptor_pdb)
                    pdbqt=df.loc[df['File name'] == filename, 'Frame'].values

                    box_center, box_min, box_max = calc_box(protein_pdb, row["residues"])
                    box_size = [abs(box_max[0] - box_min[0]), abs(box_max[1] - box_min[1]), abs(box_max[2] - box_min[2])]

                    ligands = glob.glob(os.path.join(ligand_folder, '*.pdbqt'))

                    ligand_number = len(ligands_smiles)
                    chosens = [f for f in os.listdir(prot_pdbqt_folder) if f.endswith(".pdbqt")]
                    chosen_number = len(chosens)

                    simulation_number = ligand_number * chosen_number #amount of docking simulations

                    with st.spinner("Docking simulation is now running... Please wait."):
                        for lig_path in ligands:
                            out_path = os.path.join(out_folder, os.path.basename(lig_path)[:-6] + '_smina.sdf')
                            
                            output, error = run_smina(lig_path, receptor_pdbqt, out_path, box_center, box_size, "smina")
                            st.write(output)
                            st.write(error)
                            print(output, error)
                            df_output = parse_smina_log(output)
                            df_output['library'] = ligand_folder
                            df_output['ligand'] = os.path.basename(lig_path)[:-6]
                            df_output['receptor'] = f"Frame {df['Frame'].iloc[0]}"
                            list_outputs.append(df_output)
                            time.sleep(0.5)

                            # Update progress bar
                            completed_simulations += 1
                            docking_progress_bar.progress(completed_simulations / total_simulations)

                # Combine and Save Outputs
                df_outputs = pd.concat(list_outputs, axis=0, ignore_index=True)
                df_outputs_path = os.path.join(ligand_folder, "df_outputs.csv")
                df_outputs.to_csv(df_outputs_path, index=False)
                st.success(f"Docking completed. Results saved to {df_outputs_path}")

                output_df = pd.read_csv(df_outputs_path)
                output_df = output_df.drop(columns="library")
                st.dataframe(output_df, hide_index=True, use_container_width=True)

                # Visualize the best docking pose for a ligand and receptor
                st.header("Visualize Docked Pose")

                lig_num = [f for f in os.listdir(out_folder) if f.endswith(".sdf")]
                file_count = len(lig_num)

                for i in range(file_count):
                    filename = os.path.join(out_folder, lig_num[i])
                    st.write(f"Visualizing: {filename}")

                    try:
                        with open(filename, 'r') as f:
                            docked_structure = f.read()

                            # Visualizing using Py3Dmol
                            viewer = py3Dmol.view(width=800, height=600)
                            protein_st = df['File name'].iloc[0]
                            pdb_files = os.path.join(processed_dir, "pdb_files")
                            protein_structure = os.path.join(pdb_files, f"{protein_st}.pdb")
                            st.write(protein_structure)
                            viewer = py3Dmol.view(width=800, height=600)
                            viewer.addModel(docked_structure, "sdf")
                            viewer.setStyle({'stick': {}})
                            viewer.setBackgroundColor('white')
                            viewer.zoomTo()
                            stmol.showmol(viewer, height=500, width=600)

                            viewer.addModel(protein_structure, "pdb")
                            viewer.addSurface(py3Dmol.VDW, {"opacity": "1", "color": "green"},{"hetflag": False})

                            st.write(f"Best docked pose for {lig_num[i]} visualized above.")

                    except FileNotFoundError:
                        st.warning(f"File {filename} not found.")

                else:
                    st.warning("No docked pose found.")

    else:
        st.error("Job ID not found. Please enter a valid Job ID.")