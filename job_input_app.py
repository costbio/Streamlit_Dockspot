import streamlit as st
import os
import json
from datetime import datetime
from data_analysis.config import get_config
from data_analysis.step1_extract_to_pdb import xtc_to_pdb, write_pdb_list
from data_analysis.step2_predict_pockets import merge_to_csv , run_p2rank

# Base folder for uploads
BASE_UPLOAD_FOLDER = 'streamlit_jobs'
os.makedirs(BASE_UPLOAD_FOLDER, exist_ok=True)

# Base folder for job statuses
STATUS_FOLDER = 'job_statuses'
os.makedirs(STATUS_FOLDER, exist_ok=True)

def get_status_file(job_id):
    return os.path.join(STATUS_FOLDER, f'{job_id}.json')

def load_status(job_id):
    status_file = get_status_file(job_id)
    if os.path.exists(status_file):
        with open(status_file, 'r') as f:
            return json.load(f)
    return {}

def update_status(job_id, status, step=None):
    statuses = load_status(job_id)
    statuses['status'] = status
    statuses['step'] = step
    with open(get_status_file(job_id), 'w') as f:
        json.dump(statuses, f)

# Streamlit UI
#st.title("Job Input App")

# Step 1: File Uploads
st.header("Upload Files")
xtc_file = st.file_uploader("Upload Trajectory File (*.xtc file)", type=["xtc"])
topology_file = st.file_uploader("Upload Topology File (*.gro or *.pdb file)", type=["pdb", "gro"])

# Start processing logic
def start_processing():
    if xtc_file and topology_file:
        # Generate job ID and get config paths
        job_id = datetime.now().strftime("%Y%m%d%H%M%S")
        st.session_state['job_id'] = job_id

        config = get_config(job_id)

        xtc_path = os.path.join(config['inputs_folder'], xtc_file.name)
        topology_path = os.path.join(config['inputs_folder'], topology_file.name)

        # Save uploaded files to the job-specific folder
        with open(xtc_path, "wb") as f: 
            f.write(xtc_file.getbuffer())
        with open(topology_path, "wb") as f:
            f.write(topology_file.getbuffer())

        # Initialize job status
        update_status(job_id, 'in_progress', 'Files uploaded')
        st.success(f"Your JOB ID is: {job_id}")

        # Step 1: Process XTC to PDB
        st.header("Converting XTC to PDB")
        update_status(job_id, 'in_progress', 'Converting XTC to PDB files...')
        with st.spinner("Processing XTC file..."):
            pdb_files = xtc_to_pdb(xtc_path, topology_path, config['pdb_dir'])
            filenames = [os.path.basename(path) for path in pdb_files]
            numberof_files = len(filenames)
            st.write(f"Generated {numberof_files} PDB files...")

        # Write pdb_list.ds
        update_status(job_id, 'in_progress', 'Merging PDB files...')
        pdb_list_file = write_pdb_list(config['pdb_dir'], os.path.join(config['pdb_dir'], "pdb_list.ds"))
        st.write("PDB list file created... Ready for P2Rank...")

        # Step 2: Analyze with P2Rank
        st.header("Detecting Binding Pockets")
        update_status(job_id, 'in_progress', 'Detecting binding pockets..')
        with st.spinner("Detecting binding pockets.."):
            run_p2rank(pdb_list_file, config['p2rank_processed_dir'])
            st.success("P2Rank processing completed.")

        # Merge pockets data
        pockets_csv = merge_to_csv(config['p2rank_processed_dir'], pdb_list_file)
        st.write(f"Pockets data saved to {os.path.basename(pockets_csv)}")
        update_status(job_id, 'completed')

        st.write(f"Your JOB ID is: {job_id} , please note it down.")

# Button to start the process
if st.button("Start Processing"):
    start_processing()