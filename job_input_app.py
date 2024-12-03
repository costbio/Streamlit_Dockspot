import streamlit as st
import os
from datetime import datetime
from data_analysis.step1_xtc_handling import xtc_to_pdb, write_pdb_list, run_p2rank
from data_analysis.step2_pocketscsv import merge_to_csv
from data_analysis.step3_clustering_and_repchoosing import cluster_and_save_representatives

# Create job directories
BASE_UPLOAD_FOLDER = 'streamlit_jobs'
os.makedirs(BASE_UPLOAD_FOLDER, exist_ok=True)

st.title("Job Input App")

# Step 1: File Uploads
st.header("Step 1: Upload Files")
xtc_file = st.file_uploader("Upload Trajectory File (*.xtc file)", type=["xtc"])
topology_file = st.file_uploader("Upload Topology File (*.gro or *.pdb file)", type=["pdb", "gro"])

def start_processing():
    if xtc_file and topology_file:
        # Create job folder
        job_id = datetime.now().strftime("%Y%m%d%H%M%S")
        job_folder = os.path.join(BASE_UPLOAD_FOLDER, job_id)
        os.makedirs(job_folder, exist_ok=True)

        xtc_path = os.path.join(job_folder, xtc_file.name)
        topology_path = os.path.join(job_folder, topology_file.name)

        with open(xtc_path, "wb") as f:
            f.write(xtc_file.getbuffer())
        
        with open(topology_path, "wb") as f:
            f.write(topology_file.getbuffer())
        
        st.success(f"Your JOB ID is as follows: {os.path.basename(job_folder)}")
        
        # Step 1: Process XTC to PDB
        st.header("Step 1: Converting XTC to PDB")
        with st.spinner("Processing XTC file..."):
            xtc_to_pdb_folder = os.path.join(job_folder, 'xtc_to_pdb_outputs')
            os.makedirs(xtc_to_pdb_folder, exist_ok=True)
            pdb_files = xtc_to_pdb(xtc_path, topology_path, xtc_to_pdb_folder)
            filenames = [os.path.basename(path) for path in pdb_files]
            st.write(f"Generated PDB files: {filenames}")

        # Write pdb_list.ds
        pdb_list_file = write_pdb_list(xtc_to_pdb_folder, os.path.join(xtc_to_pdb_folder, "pdb_list.ds"))
        st.write(f"PDB list file is created... Ready for P2rank...")

        # Step 2: Analyze with P2Rank
        st.header("Step 2: Detecting Binding Pockets")
        with st.spinner("Running P2Rank..."):
            run_p2rank(pdb_list_file, xtc_to_pdb_folder)
            st.success("P2Rank processing completed.")
        
        pockets_csv = merge_to_csv(xtc_to_pdb_folder, pdb_list_file)
        st.write(f"Pockets data saved to {os.path.basename(pockets_csv)}")

        st.write("Please do not forget to write this down...")
        st.write(f"Your JOB ID is as follows: {os.path.basename(job_folder)}")

# Button to start the process
if st.button("Start Processing"):
    start_processing()

