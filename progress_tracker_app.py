import streamlit as st
import os
import json

STATUS_FILE = 'job_status.json'

# Ensure the status file exists
if not os.path.exists(STATUS_FILE):
    with open(STATUS_FILE, 'w') as f:
        json.dump({}, f)

st.title("Job Progress Tracker")

# Function to load the job status
def load_status():
    with open(STATUS_FILE, 'r') as f:
        return json.load(f)

# Function to update the job status
def update_status(job_id, status, step=None):
    statuses = load_status()
    statuses[job_id] = {'status': status, 'step': step}
    with open(STATUS_FILE, 'w') as f:
        json.dump(statuses, f)

# Inquiry form
st.header("Check Job Status")
job_id = st.text_input("Enter your JOB ID:")

if st.button("Check Status"):
    statuses = load_status()
    
    if job_id in statuses:
        job_info = statuses[job_id]
        if job_info['status'] == 'in_progress':
            st.warning(f"Job `{job_id}` is still in progress. Current step: {job_info['step']}")
        elif job_info['status'] == 'completed':
            st.success(f"Job `{job_id}` has been completed!")
        else:
            st.error(f"Job `{job_id}` has an unknown status: {job_info['status']}")
    else:
        st.error(f"No job found with ID `{job_id}`.")

"""# Example function to simulate status updates (for development purposes)
if st.button("Simulate Progress Update"):
    example_job_id = "20231206153000"
    st.info(f"Simulating progress update for job `{example_job_id}`...")
    update_status(example_job_id, 'in_progress', 'Converting XTC to PDB')
    st.success("Status updated!")

if st.button("Simulate Completion"):
    example_job_id = "20231206153000"
    st.info(f"Marking job `{example_job_id}` as completed...")
    update_status(example_job_id, 'completed')
    st.success("Status updated!")"""
