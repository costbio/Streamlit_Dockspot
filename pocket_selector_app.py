import streamlit as st
import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from data_analysis.config import get_config
from data_analysis.step3_clustering_conf_creator import cluster_pockets
import nglview as nv
from IPython.display import display
import plotly.graph_objs as go



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
        pdb_location = os.path.join(processed_dir, 'pdb_files')

        if os.path.exists(csv_file_path):
            st.header("Step 2: Cluster and Generate Heatmap")
            
            clustering_depth = st.number_input("Enter clustering depth", min_value=1, max_value=10, value=3)

            if st.button("Start Clustering"):
                update_status(job_id_input, 'in_progress', 'Clustering...')
                with st.spinner("Clustering..."):
                    df_rep_pockets, full_heatmap, rep_heatmap, dataindex, repdataindex, aminoacid_list = cluster_pockets(
                        p2rank_output_folder=p2rank_dir, 
                        pdb_location=pdb_location, 
                        depth=clustering_depth,
                        plot=True
                    )
                    # Save the representative pockets to a CSV file
                    df_rep_pockets.to_csv(output_file_path, index=False)
                    st.success("Clustering completed. Representatives saved.")
                    update_status(job_id_input, 'completed', 'Clustering')

            # Step 3: Generate and Display Heatmap
            if os.path.exists(output_file_path):
                st.header("Step 3: Heatmap Visualization")

                df = pd.read_csv(output_file_path)
                
                # Full Heatmap
                if st.checkbox("Show Full Heatmap"):
                    st.write("Full dataset heatmap:")
                    reordered_matrix = full_heatmap.data2d
                    y_labels = dataindex
                    x_labels = aminoacid_list

                    fig = go.Figure(data=go.Heatmap(z=reordered_matrix, y=y_labels, x=x_labels,
                                                    hovertemplate="<b>Amino Acid:%{x}</b><br><b>%{y}</b><br>", colorscale='YlOrRd'))
                    fig.update_layout(
                        title='Interactive Heatmap with Clustering',
                        xaxis=dict(title='Amino Acids'),
                        yaxis=dict(title='Clusters'))
                    st.plotly_chart(fig)

                # Filtered Heatmap
                if st.checkbox("Show Heatmap with Representatives Only"):
                    st.write("Filtered heatmap (Representatives only):")
                    rep_reordered_matrix = rep_heatmap.data2d
                    y_labels = repdataindex
                    x_labels = aminoacid_list
                    
                    if len(x_labels) == rep_reordered_matrix.shape[1]:
                        rep_fig = go.Figure(data=go.Heatmap(z=rep_reordered_matrix, y=y_labels, x=x_labels,
                                                    hovertemplate="<b>Amino Acid:%{x}</b><br><b>%{y}</b><br>", colorscale='YlOrRd'))
                        rep_fig.update_layout(
                            title='Interactive Heatmap with Clustering',
                            xaxis=dict(title='Amino Acids'),
                            yaxis=dict(title='Clusters'))
                        st.plotly_chart(rep_fig)
                    else:
                        st.error("Mismatch in the number of residue names and columns in the heatmap matrix.")
                        a = len(aminoacid_list)
                        b = rep_reordered_matrix.shape[1]
                        st.error(f"{a},{b}")

                    
                # Step 4: Residue Viewer
                st.header("Step 4: Residue Viewer")
                test_data_name = st.text_input("Enter test data name (e.g., 'test_data_26')")

                if test_data_name:
                    # Extract the residue information for the given test data name
                    residues = df[df['test_data_name'] == test_data_name]['residues'].tolist()

                    if residues:
                        st.write(f"Residues for {test_data_name}: {residues}")
                    else:
                        st.error(f"No residues found for {test_data_name}")
                        
                    # Optionally display the protein structure
                    if st.checkbox("Show Protein Structure"):
                        try:
                            view = nv.show_structure_file(f"{pdb_location}/{test_data_name}.pdb")
                            view.add_cartoon()
                            display(view)
                        except Exception as e:
                            st.error(f"Error loading structure for {test_data_name}: {e}")

        else:
            st.error("Necessary files not found in the job folder. Please ensure the job has been processed correctly.")
    else:
        st.error(f"Job ID {job_id_input} not found. Please verify the Job ID.")
