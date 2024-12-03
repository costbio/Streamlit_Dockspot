import streamlit as st

# Define the pages
pages = {
    "Step 1: Job Input App": "job_input_app.py",
    "Step 2: Pocket Selector App": "additional_app.py",
    "Progress Tracker": "progress_tracker_app.py"
}

# Sidebar navigation
selected_page = st.sidebar.selectbox("Choose a page", list(pages.keys()))

# Import the selected page dynamically
page_file = pages[selected_page]
exec(open(page_file).read())
