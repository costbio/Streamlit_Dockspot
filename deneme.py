import pandas as pd
import py3Dmol
from stmol import showmol


# Load the CSV file
df = pd.read_csv("/home/bogrum/Streamlit_Dockspot/streamlit_jobs/20241206170251/processed/representatives.csv")

# Print out the column names
print(df.columns)

# Step 2: Iterate over each row and visualize the structures
def render_protein_structure(file_name, residues):
    # Fetch the 3D structure (using the file name to locate the PDB, for example)
    xyzview = py3Dmol.view(query=f'pdb:{file_name}')
    xyzview.setStyle({'cartoon': {'color': 'spectrum'}})
    
    # Highlight the residues
    # Convert the residues to the correct format for Py3Dmol
    resi_list = residues.split(' ')
    highlight_residues = [{'resi': resi} for resi in resi_list]
    
    # Set the style for highlighted residues
    xyzview.setStyle({'stick': {'color': 'red'}}, highlight_residues)
    
    # Show the 3D model
    showmol(xyzview, height=500, width=800)

# Step 3: Streamlit UI
import streamlit as st

st.title('Protein Structure Visualization')

# Show the protein structure for each row in the DataFrame
for _, row in df.iterrows():
    file_name = row['File name']
    residues = row['residues']
    
    # Display the file name and probability
    st.write(f"File: {file_name}, Probability: {row['probability']}")
    
    # Render the protein structure
    render_protein_structure(file_name, residues)
