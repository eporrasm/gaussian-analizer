import os
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
    
    # Extract Wiberg bond index matrix (first occurrence only)
    wiberg_start = content.find("Wiberg bond index matrix in the NAO basis:")
    if wiberg_start == -1:
        raise ValueError("No Wiberg bond index matrix found in file")
    wiberg_end = content.find("Wiberg bond index matrix in the NAO basis:", wiberg_start + 1)
    if wiberg_end == -1:  # If no second occurrence found
        wiberg_end = content.find("JK", wiberg_start)  # Or some other reliable endpoint
    wiberg_text = content[wiberg_start:wiberg_end]
    
    # Process Wiberg matrix
    num_atoms = len(npa_df)
    wiberg_matrix = np.zeros((num_atoms, num_atoms))
    current_row = -1
    cols = []  # Initialize cols list outside the loop
    
    for line in wiberg_text.split('\n'):
        if 'Atom' in line and '----' in line:
            # Column header line, get the column indices
            cols = [int(col) - 1 for col in line.split() if col.isdigit()]
            continue
            
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit() and cols:  # Check if cols is not empty
            current_row = int(parts[0]) - 1
            # Skip atom symbol and start from the numerical values
            values = []
            for val in parts[2:]:  # Start after atom number and symbol
                try:
                    values.append(float(val))
                except ValueError:
                    continue  # Skip non-numeric values
            for col, val in zip(cols, values):
                wiberg_matrix[current_row, col] = val
    
    wiberg_df = pd.DataFrame(wiberg_matrix)
    
    # Save to CSV
    npa_df.to_csv(output_dir / f'{filename}_npa.csv', index=False)
    wiberg_df.to_csv(output_dir / f'{filename}_wiberg.csv', index=False)
    
    return npa_df, wiberg_df
