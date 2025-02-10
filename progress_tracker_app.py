import streamlit as st
import os
import json

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


st.header("Check Job Status")
job_id = st.text_input("Enter your JOB ID:")

if st.button("Check Status"):
    job_info = load_status(job_id)
    
    if job_info:
        if job_info['status'] == 'in_progress':
            st.warning(f"Job `{job_id}` is still in progress. Current step: {job_info['step']}")
        elif job_info['status'] == 'completed':
            st.success(f"Job `{job_id}` has been completed!")
        else:
            st.error(f"Job `{job_id}` has an unknown status: {job_info['status']}")
    else:
        st.error(f"No job found with ID `{job_id}`.")