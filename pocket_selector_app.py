import streamlit as st
import stmol
import os
import io
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from data_analysis.config import get_config
import py3Dmol
from data_analysis.step3_clustering import cluster_pockets
import nglview as nv
from IPython.display import display
import plotly.graph_objs as go
from st_aggrid import AgGrid, GridOptionsBuilder
from datetime import datetime

# Session State Initialization

if 'clustering_done' not in st.session_state:
    st.session_state.clustering_done = False
if 'clustering_started' not in st.session_state:
    st.session_state.clustering_started = False
if 'chosen_points' not in st.session_state:
    st.session_state.chosen_points = []
if 'df_rep_pockets' not in st.session_state:
    st.session_state.df_rep_pockets = None
if 'full_heatmap' not in st.session_state:
    st.session_state.full_heatmap = None
if 'rep_heatmap' not in st.session_state:
    st.session_state.rep_heatmap = None
if 'dataindex' not in st.session_state:
    st.session_state.dataindex = None
if 'repdataindex' not in st.session_state:
    st.session_state.repdataindex = None
if 'aminoacid_list' not in st.session_state:
    st.session_state.aminoacid_list = None
if 'job_id' not in st.session_state:
    st.session_state.job_id = None
if 'config' not in st.session_state:
    st.session_state.config = None
if 'new_job_id' not in st.session_state:
    st.session_state.new_job_id = None
if 'new_config' not in st.session_state:
    st.session_state.new_config = None

# ------------------------
# Custom CSS
# ------------------------
st.markdown(
    """
    <style>
    .streamlit-expander {
        width: 700px; /* Set the width */
        max-height: 300px; /* Set the maximum height */
        overflow-y: auto; /* Add a scrollbar if content overflows */
    }
            
    .fixed-column {
        height: 600px;  /* Set a fixed height for the second column */
        overflow-y: auto;  /* Make it scrollable if content exceeds the height */
    }
    </style>
    """, 
    unsafe_allow_html=True
)

# ------------------------
# Base folders and helper functions
# ------------------------
BASE_JOB_FOLDER = 'streamlit_jobs'
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



def render_structure_with_residues(structure_file_path, residues_list):
    # Expand the ~ to the full path if needed
    structure_file_path = os.path.expanduser(structure_file_path)
    # Read the content of the PDB file
    with open(structure_file_path, 'r') as file:
        pdb_data = file.read()
    # Create the viewer using the PDB data as a string
    viewer = py3Dmol.view(width=700, height=500)
    viewer.addModel(pdb_data, "pdb")
    # Set the entire structure to cartoon style
    viewer.setStyle({'cartoon': {'color': 'spectrum'}})
    # Add a surface for the selected residues (e.g. in red)
    viewer.addSurface(py3Dmol.VDW, {"opacity": 0.8, "color": "red"}, {"serial": residues_list})
    # Zoom to the structure
    viewer.zoomTo()
    # Show the viewer using stmol
    stmol.showmol(viewer, height=500, width=600)
    
    # Format residues list in a table
    items_per_row = 6
    rows = [residues_list[i:i + items_per_row] for i in range(0, len(residues_list), items_per_row)]
    columns = [f"Atom {i + 1}" for i in range(items_per_row)]
    df = pd.DataFrame(rows, columns=columns)
    df.index = [""] * len(df)
    with st.expander("View Atom Numbers"):
        st.table(df)

# ------------------------
# UI: Job ID Input
# ------------------------
st.header("Enter Your Job ID")
job_id_input = st.text_input("Enter your Job ID", value=st.session_state.get("job_id", ""))

if job_id_input:
    # Save the user's input as the "good old" job ID
    st.session_state.job_id = job_id_input
    job_folder = os.path.join(BASE_JOB_FOLDER, st.session_state.job_id)
    if os.path.exists(job_folder):
        st.success(f"Job found: {job_id_input}")
        
        # Generate a NEW job id and config for outputs if not already created
        if st.session_state.new_job_id is None:
            new_job_id = datetime.now().strftime("%Y%m%d%H%M%S")
            st.session_state.new_job_id = new_job_id
            st.session_state.new_config = get_config(new_job_id)
        st.write(f"Your new Job ID is: {st.session_state.new_job_id}")
        
        # Define folder paths using the input (old) job id
        processed_dir = os.path.join(job_folder, 'processed')
        p2rank_dir = os.path.join(processed_dir, "p2rank_output")
        csv_file_path = os.path.join(p2rank_dir, 'pockets.csv')
        # Output file uses the new job id's configuration
        output_file_path = os.path.join(st.session_state.new_config['processed_dir'], 'representatives.csv')
        pdb_location = os.path.join(processed_dir, 'pdb_files')

        if os.path.exists(csv_file_path):
            st.header("Cluster and Generate Heatmap")
            clustering_depth = st.number_input("Enter clustering depth", min_value=1, max_value=10, value=3)

            # ------------------------
            # Clustering Button
            # ------------------------
            if st.button("Start Clustering"):
                st.session_state.clustering_started = True
                update_status(st.session_state.job_id, 'in_progress', 'Clustering...')
                with st.spinner("Clustering..."):
                    df_rep_pockets, full_heatmap, rep_heatmap, dataindex, repdataindex, aminoacid_list, surface_atoms = cluster_pockets(
                        p2rank_output_folder=p2rank_dir, 
                        pdb_location=pdb_location, 
                        depth=clustering_depth,
                        plot=True
                    )

                    # Save results into session state
                    st.session_state.df_rep_pockets = df_rep_pockets
                    st.session_state.full_heatmap = full_heatmap
                    st.session_state.rep_heatmap = rep_heatmap
                    st.session_state.dataindex = dataindex
                    st.session_state.repdataindex = repdataindex
                    st.session_state.aminoacid_list = aminoacid_list
                    st.session_state.surface_atoms = surface_atoms  # Will be updated later

                    st.session_state.clustering_done = True

                    # Save representatives CSV using new job id configuration
                    df_rep_pockets.to_csv(output_file_path, index=False)
                    st.success("Clustering completed. Representatives saved.")
                    update_status(st.session_state.job_id, 'completed', 'Clustering')

            # ------------------------
            # Display Results if Clustering Completed
            # ------------------------
            if st.session_state.clustering_done and os.path.exists(output_file_path):
                st.header("Found Pockets")
                df = pd.read_csv(output_file_path)
                
                # Full Heatmap
                if st.checkbox("Show Full Heatmap"):
                    st.write("Full dataset heatmap:")
                    reordered_matrix = st.session_state.full_heatmap.data2d
                    y_labels = st.session_state.dataindex
                    x_labels = st.session_state.aminoacid_list

                    fig = go.Figure(
                        data=go.Heatmap(
                            z=reordered_matrix,
                            y=y_labels,
                            x=x_labels,
                            hovertemplate="<b>Amino Acid: %{x}</b><br><b>%{y}</b><br>",
                            colorscale='YlOrRd'
                        )
                    )
                    fig.update_layout(
                        title='Interactive Heatmap with Clustering',
                        xaxis=dict(title='Amino Acids'),
                        yaxis=dict(title='Clusters')
                    )
                    st.plotly_chart(fig)

                # Two columns: one for heatmap, one for representative selection
                col1, col2 = st.columns([3, 1])

                # Left Column: Filtered (Representative) Heatmap
                with col1:
                    st.write("Filtered heatmap (Representatives only):")
                    rep_reordered_matrix = st.session_state.rep_heatmap.data2d
                    y_labels = st.session_state.repdataindex
                    x_labels = st.session_state.aminoacid_list
                    
                    if len(x_labels) == rep_reordered_matrix.shape[1]:
                        rep_fig = go.Figure(
                            data=go.Heatmap(
                                z=rep_reordered_matrix,
                                y=y_labels,
                                x=x_labels,
                                hovertemplate="<b>Amino Acid: %{x}</b><br><b>%{y}</b><br>",
                                colorscale='YlOrRd'
                            )
                        )
                        rep_fig.update_layout(
                            title='Interactive Heatmap with Clustering',
                            xaxis=dict(title='Amino Acids'),
                            yaxis=dict(title='Clusters')
                        )
                        st.plotly_chart(rep_fig)
                    else:
                        st.error("Mismatch in the number of residue names and columns in the heatmap matrix.")
                        st.error(f"{len(x_labels)} residue names vs {rep_reordered_matrix.shape[1]} columns.")

                # Right Column: Representative Selection Form
                with col2:
                    st.write("Choose a conformation:")
                    # Use a form to capture checkbox selections together
                    with st.form("selection_form"):
                        form_chosen_points = []
                        # Iterate over the representative labels in reverse order
                        for label in st.session_state.get('repdataindex', [])[::-1]:
                            if st.checkbox(f"Select {label}", key=f"checkbox_{label}"):
                                form_chosen_points.append(label)
                        submitted = st.form_submit_button("Save selections")
                        if submitted:
                            st.session_state.chosen_points = form_chosen_points
                            st.success("Selections saved!")
                            
                
                # Display structures if any selections have been made
                if st.session_state.chosen_points:
                    st.header("Structures for Selected Representatives")
                    chosen_points = st.session_state.chosen_points
                    num_chosen_points = len(chosen_points)
                    # Loop in pairs so you get two columns per row
                    for i in range(0, num_chosen_points, 2):
                        col_left, col_right = st.columns(2)
                        # Left structure
                        with col_left:
                            label_1 = chosen_points[i]
                            #atoms = df.loc["Frame:35 Index:1", "surf_atom_ids"]

                            df_rep = st.session_state.df_rep_pockets
                            # Retrieve the filename corresponding to the label
                            label_1_filename = df_rep.loc[df_rep.index.isin([label_1]), 'File name'].iloc[0]
                            pdb_file_1 = os.path.join(pdb_location, f"{label_1_filename}.pdb")
                            if os.path.exists(pdb_file_1):
                                st.subheader(f"{label_1}")
                                # Assume residue numbers are extracted from the column names (modify as needed)
                                nums = [col.split("_")[1] for col in df_rep.drop(
                                    columns=['File name', 'Frame', 'pocket_index', 'probability', 'residues', "frame_pocket"]
                                ).columns]
                                residues = st.session_state.aminoacid_list
                                atom_numbers = st.session_state.surface_atoms[chosen_points[i]]

                                render_structure_with_residues(pdb_file_1, atom_numbers)
                            else:
                                st.warning(f"PDB file for {label_1} not found. Path was {pdb_file_1}")
                        # Right structure (if exists)
                        if i + 1 < num_chosen_points:
                            with col_right:
                                label_2 = chosen_points[i + 1]
                                
                                df_rep = st.session_state.df_rep_pockets
                                label_2_filename = df_rep.loc[df_rep.index.isin([label_2]), 'File name'].iloc[0]
                                pdb_file_2 = os.path.join(pdb_location, f"{label_2_filename}.pdb")
                                if os.path.exists(pdb_file_2):
                                    st.subheader(f"{label_2}")
                                    nums = [col.split("_")[1] for col in df_rep.drop(
                                        columns=['File name', 'Frame', 'pocket_index', 'probability', 'residues', "frame_pocket"]
                                    ).columns]
                                    atom_numbers2 = st.session_state.surface_atoms[chosen_points[i]]
                                    render_structure_with_residues(pdb_file_2,atom_numbers2)
                                else:
                                    st.warning(f"PDB file for {label_2} not found. Path was {pdb_file_2}")

                # Button to save the selected representatives to CSV
                if st.button("Save selections"):
                    if st.session_state.chosen_points:
                        csv_path = os.path.join(st.session_state.new_config['processed_dir'], "chosen_representatives.csv")
                        df_rep = st.session_state.df_rep_pockets
                        # Filter the dataframe for the selected labels
                        filtered_df = df_rep.loc[st.session_state.chosen_points]
                        filtered_df.to_csv(csv_path, index=False)
                        st.success("Selections are saved to a CSV file.")
                    else:
                        st.warning("No representatives selected.")
        else:
            st.error("Necessary files not found in the job folder. Please ensure the job has been processed correctly.")
    else:
        st.error(f"Job ID {job_id_input} not found. Please verify the Job ID.")
