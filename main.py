import streamlit as st
from streamlit_option_menu import option_menu
import extra_streamlit_components as stx
import os
from streamlit_extras.app_logo import add_logo

st.set_page_config(
    page_title="DockSpot",
    page_icon="logo_2.png",  # You can add an emoji or path to an icon
    layout="wide",  # This spreads the app across the screen
    initial_sidebar_state="collapsed",  # Collapse sidebar initially
    menu_items={
        'Get Help': 'https://www.extremelycoolapp.com/help',
        'Report a bug': "https://www.extremelycoolapp.com/bug",
        'About': "# This is a header. This is an *extremely* cool app!"
    }
)
#st.logo("logo_2.png",size="large")
#add_logo("/home/bogrum/Streamlit_Dockspot/logo.png", height=300)
st.image('logo.png', width=200)

# Define the pages
pages = {
    "PDB Generator": "job_input_app.py",
    "Pocket Finder": "pocket_selector_app.py",
    "Docking Simulator": "docking.py",
    "Task Monitor": "progress_tracker_app.py"
}

# Horizontal menu
selected = option_menu(
    None, 
    list(pages.keys()), 
    icons=['box-arrow-in-right', 'hand-index', "gear", 'percent'], 
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
