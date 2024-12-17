import streamlit as st
from streamlit_option_menu import option_menu
import os

# Define the pages
pages = {
    "Job Input App": "job_input_app.py",
    "Pocket Selector App": "pocket_selector_app.py",
    "Docking App": "docking.py",
    "Progress Tracker": "progress_tracker_app.py"
}

# Horizontal menu
selected = option_menu(
    None, 
    list(pages.keys()), 
    icons=['house', 'cloud-upload', "list-task", 'gear'], 
    menu_icon="cast", 
    default_index=0, 
    orientation="horizontal"
)

# Get the corresponding file
selected_file = pages[selected]

# Execute the selected file
with open(selected_file) as f:
    code = f.read()
    exec(code)
