import pandas as pd
import numpy as np
from pathlib import Path

def read_gaussian_matrices(filename):
    """
    Reads a Gaussian log file and extracts the Natural Population Analysis and 
    Wiberg bond index matrices.
    
    Args:
        filename (str): Name of the log file in the input folder
    
    Returns:
        tuple: (npa_df, wiberg_df) - Two pandas DataFrames containing the matrices
    """
    return read_npa_matrix(filename), read_wiberg_matrix(filename)

def read_npa_matrix(filename):
    """
    Reads a Natural Population Analysis table from a Gaussian log file.
    
    Args:
        filename (str): Name of the log file in the input folder
    
    Returns:
        pd.DataFrame: DataFrame containing the NPA table
    """
    # Define input/output paths
    input_dir = Path('input')
    output_dir = Path('output')
    file_path = input_dir / filename
    
    # Read the file
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Extract Natural Population Analysis (first occurrence only)
    npa_start = content.find("Summary of Natural Population Analysis:")
    if npa_start == -1:
        raise ValueError("No Natural Population Analysis found in file")
    npa_end = content.find("=======================================================================", npa_start)
    npa_text = content[npa_start:npa_end]

    # Process NPA into DataFrame
    # Find the actual data start after the headers
    data_start = npa_text.find(" -------------------------------------------")
    npa_data = npa_text[data_start:].split('\n')
    # Skip header lines and empty lines, only take lines with actual data
    npa_lines = []
    header_found = False
    for line in npa_data:
        if 'Atom  No' in line:  # Skip the header line
            header_found = True
            continue
        if header_found and line.strip() and not line.startswith(' ---'):
            parts = line.split()
            if len(parts) >= 7:
                try:
                    npa_lines.append([
                        parts[0],  # Keep full atom label
                        parts[1],  # Atom number
                        float(parts[2]),  # Natural Charge
                        float(parts[3]),  # Core
                        float(parts[4]),  # Valence
                        float(parts[5]),  # Rydberg
                        float(parts[6])   # Total
                    ])
                except ValueError:
                    continue  # Skip lines that can't be properly converted
    
    npa_df = pd.DataFrame(npa_lines, 
                         columns=['Atom_Label', 'Atom_Num', 'Natural_Charge', 'Core', 'Valence', 'Rydberg', 'Total'])

    # Convert only numeric columns to float
    numeric_columns = ['Natural_Charge', 'Core', 'Valence', 'Rydberg', 'Total']
    for col in numeric_columns:
        npa_df[col] = pd.to_numeric(npa_df[col], errors='coerce')
    
    # Keep Atom_Label and Atom_Num as strings
    npa_df['Atom_Label'] = npa_df['Atom_Label'].astype(str)
    npa_df['Atom_Num'] = npa_df['Atom_Num'].astype(str)

    # Save to output folder
    npa_df.to_csv(output_dir / f'{filename}_npa.csv', index=False)

    return npa_df

def read_wiberg_matrix(filename):
    return None