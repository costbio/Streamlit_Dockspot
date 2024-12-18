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
from data_analysis.step3_clustering_conf_creator import cluster_pockets
import nglview as nv
from IPython.display import display
import plotly.graph_objs as go
from st_aggrid import AgGrid, GridOptionsBuilder

st.markdown("""
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
""", unsafe_allow_html=True)

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

def render_structure_with_residues(structure_file_path, res_list, residues_list):
    import py3Dmol
    import stmol

    # Expand the ~ to the full path if needed (this is optional depending on your input)
    import os
    structure_file_path = os.path.expanduser(structure_file_path)

    # Read the content of the PDB file
    with open(structure_file_path, 'r') as file:
        pdb_data = file.read()

    # Create the viewer using the PDB data as a string
    viewer = py3Dmol.view(width=700, height=500)
    viewer.addModel(pdb_data, "pdb")
    viewer.addSurface(py3Dmol.VDW,{"opacity": 0.5, "color": "blue"})  # Set the default style to cartoon

    viewer.addSurface(py3Dmol.VDW,{"opacity": 0.8, "color": "red"},{"resi": res_list})


    # Zoom to fit the structure
    viewer.zoomTo()

    # Show the viewer using stmol
    stmol.showmol(viewer, height=500, width=600)
    
    
    # Define the number of items per row
    items_per_row = 6

    # Format the data into rows
    rows = [residues_list[i:i + items_per_row] for i in range(0, len(residues_list), items_per_row)]

    # Create a DataFrame with meaningful column headers
    columns = [f"Residue {i + 1}" for i in range(items_per_row)]
    df = pd.DataFrame(rows, columns=columns)

    # Hide the index in the table by resetting it to blank strings
    df.index = [""] * len(df)

    # Display the table
    with st.expander("View Residues"):
        st.table(df)

# Streamlit UI
#st.title("Pocket Selector App")

# Step 1: Job ID Input
st.header("Enter Your Job ID")
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
            st.header("Cluster and Generate Heatmap")
            
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

                    # Save clustering results to session state
                    st.session_state['df_rep_pockets'] = df_rep_pockets
                    st.session_state['full_heatmap'] = full_heatmap
                    st.session_state['rep_heatmap'] = rep_heatmap
                    st.session_state['dataindex'] = dataindex
                    st.session_state['repdataindex'] = repdataindex
                    st.session_state['aminoacid_list'] = aminoacid_list
                    st.session_state['clustering_done'] = True

                    # Save the representative pockets to a CSV file
                    df_rep_pockets.to_csv(output_file_path, index=False)
                    st.success("Clustering completed. Representatives saved.")
                    update_status(job_id_input, 'completed', 'Clustering')

            # Step 3: Generate and Display Heatmap
            if os.path.exists(output_file_path):
                st.header("Found Pockets")

                df = pd.read_csv(output_file_path)
                
                # Full Heatmap
                if st.checkbox("Show Full Heatmap"):
                    st.write("Full dataset heatmap:")
                    reordered_matrix = st.session_state['full_heatmap'].data2d
                    y_labels = st.session_state['dataindex']
                    x_labels = st.session_state['aminoacid_list']

                    fig = go.Figure(data=go.Heatmap(z=reordered_matrix, y=y_labels, x=x_labels,
                                                    hovertemplate="<b>Amino Acid:%{x}</b><br><b>%{y}</b><br>", colorscale='YlOrRd'))
                    fig.update_layout(
                        title='Interactive Heatmap with Clustering',
                        xaxis=dict(title='Amino Acids'),
                        yaxis=dict(title='Clusters'))
                    st.plotly_chart(fig)

                # Create two columns
                col1, col2 = st.columns([3, 1])

                with col1:
                    # Filtered heatmap
                    st.write("Filtered heatmap (Representatives only):")
                    rep_reordered_matrix = st.session_state['rep_heatmap'].data2d
                    y_labels = st.session_state['repdataindex']
                    x_labels = st.session_state['aminoacid_list']
                    chosen_points = []
                    
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
                        a = len(x_labels)
                        b = rep_reordered_matrix.shape[1]
                        st.error(f"{a},{b}")

                # Step 4: Heatmap rep chooser.
                # Display a table for representatives with checkboxes
                with col2:
                    st.write("Choose a conformation:")
                    with st.container():
                        # Apply the scrollable class to make the column content scrollable
                        st.markdown('<div class="stScrollBox">', unsafe_allow_html=True)

                        # Create checkboxes for each row label
                        chosen_points = []
                        for i, label in enumerate(y_labels[::-1]):
                            checkbox = st.checkbox(f"Select {label}", key=f"checkbox_{i}")
                            if checkbox:
                                chosen_points.append(label)

                        st.markdown('</div>', unsafe_allow_html=True)  # Closing div tag=True

                if chosen_points:
                    st.header("Structures for Selected Representatives")

                    # Create two columns for each pair of structures
                    num_chosen_points = len(chosen_points)
                    for i in range(0, num_chosen_points, 2):
                        # Ensure there's no index out of range error when creating columns
                        col1, col2 = st.columns(2)
                        
                        # For the first column
                        with col1:
                            label_1 = chosen_points[i]
                            df_rep_pockets = st.session_state['df_rep_pockets']
                            df = df_rep_pockets
                            label_1_filename = df.loc[df.index.isin([label_1]), 'File name'].iloc[0]
                            pdb_file_1 = os.path.join(pdb_location, f"{label_1_filename}.pdb")
                            
                            if os.path.exists(pdb_file_1):
                                st.subheader(f"{label_1}")
                                nums = [col.split("_")[1] for col in df_rep_pockets.drop(columns=['File name','Frame','pocket_index','probability','residues',"frame_pocket"]).columns]
                                residues=st.session_state['aminoacid_list']
                                render_structure_with_residues(pdb_file_1, nums,residues)
                            else:
                                st.warning(f"PDB file for {label_1} not found. Path was {pdb_file_1}")

                        # For the second column, only if there's a second item
                        if i + 1 < num_chosen_points:
                            with col2:
                                label_2 = chosen_points[i + 1]
                                label_2_filename = df.loc[df.index.isin([label_2]), 'File name'].iloc[0]
                                pdb_file_2 = os.path.join(pdb_location, f"{label_2_filename}.pdb")
                                
                                if os.path.exists(pdb_file_2):
                                    st.subheader(f"{label_2}")
                                    nums = [col.split("_")[1] for col in df_rep_pockets.drop(columns=['File name','Frame','pocket_index','probability','residues',"frame_pocket"]).columns]
                                    render_structure_with_residues(pdb_file_2, nums,st.session_state['aminoacid_list'])
                                else:
                                    st.warning(f"PDB file for {label_2} not found. Path was {pdb_file_2}")

                
                

                    # Save selected checkboxes to a text file
                    if st.button("Save selections"):
                        if chosen_points:
                            csv_path = os.path.join(processed_dir, "chosen_representatives.csv")
                            with open(csv_path, "w") as f:
                                for point in chosen_points:
                                    df_rep_pockets = st.session_state['df_rep_pockets']
                                    filtered_df = df_rep_pockets.loc[chosen_points]

                                    # Save the filtered DataFrame to a CSV file
                                    #csv_path = os.path.join(job_folder, "chosen_representatives.csv")
                                    filtered_df.to_csv(csv_path, index=False)

                                    st.success("Selections are saved to a dataframe.")
                        else:
                            st.warning("No representatives selected.")
        else:
            st.error("Necessary files not found in the job folder. Please ensure the job has been processed correctly.")
    else:
        st.error(f"Job ID {job_id_input} not found. Please verify the Job ID.")
