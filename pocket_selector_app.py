import streamlit as st
import os
import json
from datetime import datetime
from data_analysis.config import get_config
from data_analysis.step3_clustering_and_repchoosing import cluster_and_save_representatives
import pandas as pd
from data_analysis.step4_proteinviewer import load_and_process_data, get_residues_for_test_data

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

# Streamlit UI
st.title("Pocket Selector App")

# Step 1: Job ID Input
st.header("Step 1: Enter Your Job ID")
job_id_input = st.text_input("Enter your Job ID")

# Check if the job exists
if job_id_input:
    job_folder = os.path.join(BASE_JOB_FOLDER, job_id_input)
    if os.path.exists(job_folder):
        st.success(f"Job found: {job_id_input}")
        
        # Get paths from the job folder
        processed_dir = os.path.join(job_folder, 'processed')
        p2rank_dir = os.path.join(processed_dir, "p2rank_output")
        csv_file_path = os.path.join(p2rank_dir, 'pockets.csv')
        output_file_path = os.path.join(processed_dir, 'representatives.csv')

        # Check if the necessary files exist
        if os.path.exists(csv_file_path):
            st.write("Found the necessary files for clustering analysis.")
            
            # Step 2: Clustering and Save Representatives
            st.header("Step 2: Cluster and Save Representative Pockets")
            clustering_depth = st.number_input("Enter clustering depth", min_value=1, max_value=10, value=3)

            if st.button("Start Clustering"):
                update_status(job_id_input, 'in_progress', 'Clustering...')
                with st.spinner("Clustering..."):
                    representatives = cluster_and_save_representatives(
                        csv_file_path, output_file_path, clustering_depth
                    )
                    st.success("Clustering completed. Representatives saved.")
                    st.write(f"Representative pockets saved to: {output_file_path}")
                    update_status(job_id_input, 'completed')

            # Step 3: Load and Process Data (Step 4 Integration)
            st.header("Step 3: View Residues for Test Data")

            if os.path.exists(output_file_path):
                st.write("Loading processed data...")
                df = load_and_process_data(output_file_path)

                # Allow user to input a test data name
                test_data_name = st.text_input("Enter test data name (e.g., 'test_data_26')")

                if test_data_name:
                    residues = get_residues_for_test_data(df, test_data_name)

                    if residues is not None:
                        st.write(f"Residues for {test_data_name}: {residues}")
                    else:
                        st.error(f"No residues found for {test_data_name}")
                
                # Display the first few rows of the DataFrame
                st.write(df)

        else:
            st.error("Necessary files not found in the job folder. Please ensure the job has been processed correctly.")
    else:
        st.error(f"Job ID {job_id_input} not found. Please verify the Job ID.")
