import streamlit as st
import os
import json
import pandas as pd
from data_analysis.step4_docking import *
from data_analysis.step4_ligand_handling import smiles_to_pdbqt
import time
import py3Dmol
import numpy as np

# Base folder for jobs
BASE_JOB_FOLDER = 'streamlit_jobs'

# Job status file
STATUS_FILE = 'job_status.json'

# Ensure the status file exists
if not os.path.exists(STATUS_FILE):
    with open(STATUS_FILE, 'w') as f:
        json.dump({}, f)

# Functions for job tracking
def load_status():
    with open(STATUS_FILE, 'r') as f:
        return json.load(f)

def update_status(job_id, status, step=None):
    statuses = load_status()
    statuses[job_id] = {'status': status, 'step': step}
    with open(STATUS_FILE, 'w') as f:
        json.dump(statuses, f)

# Function to clean the pdbqt_files folder
def clean_pdbqt_folder(folder_path):
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if os.path.isfile(file_path):
            os.remove(file_path)

# Streamlit UI
#st.title("Docking App")

# Step 1: Job ID Input
st.header("Enter Your Job ID")
job_id_input = st.text_input("Enter your Job ID")

# Check if the job exists
if job_id_input:
    job_folder = os.path.join(BASE_JOB_FOLDER, job_id_input)
    if os.path.exists(job_folder):
        st.success(f"Job found: {job_id_input}")

        # Display job status
        job_status = load_status().get(job_id_input, {}).get('status', 'Not Started')
        #st.write(f"Job Status: {job_status}")

        # Step 2: Ligand SMILES Input
        st.header("Enter SMILES")
        smiles_input = st.text_area("Enter SMILES strings for ligands: (Enter one SMILES string per line, press Enter after each one):", 
                                    help="Enter SMILES strings for the ligands, separated by lines.")
        
        # Process Ligands Button
        if st.button("Process Ligands"):
            if not smiles_input:
                st.warning("Please enter SMILES strings to proceed.")
            else:
                ligands_smiles = [smiles.strip() for smiles in smiles_input.split('\n')]
                df = pd.DataFrame(ligands_smiles, columns=["Ligands"])
                st.dataframe(df, hide_index=True,use_container_width=True)
                st.info("Converting ligands...")

                # Creating necessary folders
                processed_dir = os.path.join(job_folder, "processed")
                pdbqt_folder = os.path.join(processed_dir, "pdbqt_files")
                os.makedirs(pdbqt_folder, exist_ok=True)

                # Clean previous pdbqt files
                clean_pdbqt_folder(pdbqt_folder)
                ligand_folder = os.path.join(pdbqt_folder, "ligand_pdbqt_files")
                os.makedirs(ligand_folder, exist_ok=True)

                prot_pdbqt_folder = os.path.join(pdbqt_folder,"protein_pdbqt_files")
                os.makedirs(prot_pdbqt_folder, exist_ok=True)

                # Loop through SMILES and process each ligand
                for i, smiles in enumerate(ligands_smiles):
                    pdbqt_path = os.path.join(ligand_folder, f"ligand_{i + 1}.pdbqt")
                    result_message = smiles_to_pdbqt(smiles, pdbqt_path)
                    #st.write(f"Ligand {i + 1} Conversion Status: Converted successfully.")
                    time.sleep(1)  # simulate processing time

                st.success("All ligands are processed successfully!")


                # Step 3: Docking Setup
                st.header("Ensemble Docking")
                # Load the chosen representatives dataframe
                csv_path = os.path.join(processed_dir, "chosen_representatives.csv")
                df = pd.read_csv(csv_path)
                new_df = df[["Frame","pocket_index","probability","residues"]]
                st.write("Chosen receptors")
                st.dataframe(new_df,hide_index=True, use_container_width=True)
                
                out_folder = os.path.join(processed_dir, "docking_results")
                os.makedirs(out_folder, exist_ok=True)

                # Iterate through receptors and ligands
                list_outputs = []
                for index, row in df.iterrows():
                    receptor_pdb = os.path.join(processed_dir, "pdb_files", f"{row['File name']}.pdb")
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
                    #print(f"Looking for file: {filename}")
                    #print(df['File name'].head())

                    #st.write(f"Receptor PDBQT file is created successfully for: {pdbqt[0]}")

                    

                    box_center, box_min, box_max = calc_box(protein_pdb, row["residues"])
                    box_size = [abs(box_max[0] - box_min[0]), abs(box_max[1] - box_min[1]), abs(box_max[2] - box_min[2])]

                    ligands = glob.glob(os.path.join(ligand_folder, '*.pdbqt'))

                    ligand_number = len(ligands_smiles)
                    chosens = [f for f in os.listdir(prot_pdbqt_folder) if f.endswith(".pdbqt")]
                    chosen_number = len(chosens)


                    simulation_number = ligand_number * chosen_number #amount of docking simulations
                    #st.write(simulation_number)

                    #if simulation_number > 2:
                    #   st.write("Simulation number is within limits. Please do run our customized script to carry on with your analysis.")
                    #  st.stop()


                    for lig_path in ligands:
                        out_path = os.path.join(out_folder, os.path.basename(lig_path)[:-6] + '_smina.sdf')
                        with st.spinner("Docking simulation is now running... Please wait."):
                            output = run_smina(lig_path, receptor_pdbqt, out_path, box_center, box_size, "smina")
                        with st.spinner("Docking affinities are being saved... Please wait."):
                            df_output = parse_smina_log(output)
                            df_output['library'] = ligand_folder
                            df_output['ligand'] = os.path.basename(lig_path)[:-6]
                            df_output['receptor'] = os.path.basename(receptor_pdb)
                            list_outputs.append(df_output)

                # Combine and Save Outputs
                df_outputs = pd.concat(list_outputs, axis=0, ignore_index=True)
                df_outputs_path = os.path.join(ligand_folder, "df_outputs.csv")
                df_outputs.to_csv(df_outputs_path, index=False)
                st.success(f"Docking completed. Results saved to {df_outputs_path}")

                # Visualize the best docking pose for a ligand and receptor
                st.header("Visualize Docked Pose")

                lig_num = [f for f in os.listdir(out_folder) if f.endswith(".sdf")]
                file_count = len(lig_num)

                for i in range(file_count):
                    filename= os.path.join(out_folder,lig_num[i])
                    st.write(filename)
                    with open(filename, 'r') as f:
                        docked_structure = f.read()

                        # Visualizing using Py3Dmol
                        viewer = py3Dmol.view(width=800, height=600)
                        viewer.addModel(docked_structure, "sdf")
                        viewer.setStyle({'stick': {}})
                        viewer.setBackgroundColor('white')
                        viewer.zoomTo()
                        viewer.show()
                        st.write("Best docked pose visualized above.")

                else:
                    st.warning("No docked pose found.")

    else:
        st.error("Job ID not found. Please enter a valid Job ID.")
