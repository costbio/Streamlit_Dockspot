import pandas as pd

def read_representatives(file_path):
    data = []
    current_data = {}
    
    # Open and read the CSV file
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            
            if not line:
                continue
            
            parts = line.split(",", 1)
            if len(parts) == 2:
                key, value = parts
                current_data[key] = value
            
            if 'residues' in current_data:
                data.append(current_data)
                current_data = {}

    return data

def process_residues_column(df):
    df['residues'] = df['residues'].apply(lambda x: x.split() if isinstance(x, str) else [])
    return df

def load_and_process_data(file_path):
    data = read_representatives(file_path)
    df = pd.DataFrame(data)
    df = process_residues_column(df)
    return df

def get_residues_for_test_data(df, test_data_name):
    """
    Filters the DataFrame to get the residues for a specific 'File name'.
    
    Args:
        df (pd.DataFrame): The DataFrame with processed data.
        test_data_name (str): The name of the test data to filter by.
        
    Returns:
        list: The residues associated with the given 'File name', or None if no match is found.
    """
    result = df[df['File name'] == test_data_name]
    if not result.empty:
        return result['residues'].iloc[0]  # Return the first match's residues
    else:
        return None  # Return None if no match is found
