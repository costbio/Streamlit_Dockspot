import os
import glob

def get_config(job_id):
    """
    Returns configuration paths for the given job ID.
    """
    # Base job folder
    BASE_UPLOAD_FOLDER = os.path.join(os.getcwd(), 'streamlit_jobs')

    # Paths specific to the job ID
    job_folder = os.path.join(BASE_UPLOAD_FOLDER, job_id)
    inputs_folder = os.path.join(job_folder, 'inputs')
    processed_dir = os.path.join(job_folder, 'processed')

    # Ensure directories exist
    os.makedirs(inputs_folder, exist_ok=True)
    os.makedirs(processed_dir, exist_ok=True)

    # File paths for PDB and P2Rank outputs
    pdb_dir = os.path.join(processed_dir, 'pdb_files')
    os.makedirs(pdb_dir, exist_ok=True)

    p2rank_processed_dir = os.path.join(processed_dir, 'p2rank_output')
    os.makedirs(p2rank_processed_dir, exist_ok=True)

    return {
        'job_folder': job_folder,
        'inputs_folder': inputs_folder,
        'processed_dir': processed_dir,
        'pdb_dir': pdb_dir,
        'p2rank_processed_dir': p2rank_processed_dir
    }

# Optional: Function to retrieve file paths dynamically based on job folder
def get_input_files(inputs_folder):
    """
    Reads input file paths from a job's input folder.
    Assumes text files in the inputs folder contain specific file paths.
    """
    input_files = {}
    for file_name in ['xtc_directory.txt', 'sdf_directory.txt', 'pdb_directory.txt', 'topology_directory.txt']:
        file_path = os.path.join(inputs_folder, file_name)
        if os.path.exists(file_path):
            with open(file_path, 'r') as f:
                input_files[file_name] = f.read().strip()
    return input_files
