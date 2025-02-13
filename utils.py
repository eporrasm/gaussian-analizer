import pandas as pd
import numpy as np
from pathlib import Path
from config import INPUT_DIR, OUTPUT_DIR, OUTPUT_EXTENSION, NPA_START_MARKER, NPA_END_MARKER, WIBERG_START_MARKER, WIBERG_END_MARKER, NPA_COLUMNS

def read_gaussian_matrices(filename):
    """
    Reads a Gaussian log file and extracts the Natural Population Analysis and 
    Wiberg bond index matrices.
    
    Args:
        filename (str): Name of the log file in the input folder
    
    Returns:
        tuple: (pd.DataFrame, np.ndarray) - Pandas DataFrame containing the NPA matrix
               and Numpy array containing the Wiberg bond index matrix
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
    input_dir = Path(INPUT_DIR)
    output_dir = Path(OUTPUT_DIR)
    file_path = input_dir / filename
    
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Extract Natural Population Analysis
    npa_start = content.find(NPA_START_MARKER)
    if npa_start == -1:
        raise ValueError("No Natural Population Analysis found in file")
    npa_end = content.find(NPA_END_MARKER, npa_start)
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
    
    df = pd.DataFrame(npa_lines, columns=NPA_COLUMNS)

    # Convert only numeric columns to float
    numeric_columns = ['Natural_Charge', 'Core', 'Valence', 'Rydberg', 'Total']
    for col in numeric_columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Keep Atom_Label and Atom_Num as strings
    df['Atom_Label'] = df['Atom_Label'].astype(str)
    df['Atom_Num'] = df['Atom_Num'].astype(str)

    # Save to output folder
    output_path = output_dir / f'{filename}_npa{OUTPUT_EXTENSION}'
    df.to_csv(output_path, index=False)

    return df

def read_wiberg_matrix(filename):
    """
    Reads a Wiberg bond index matrix from a Gaussian log file.
    The matrix is split into multiple blocks.
    
    Args:
        filename (str): Name of the log file in the input folder
    
    Returns:
        np.ndarray: Numpy array containing the Wiberg bond index matrix
    """
    input_dir = Path(INPUT_DIR)
    output_dir = Path(OUTPUT_DIR)
    file_path = input_dir / filename

    matrix = None
    found_wiberg = False
    current_cols = None
    reading_data = False
    n = get_wiberg_matrix_size(filename)
    
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            
            if WIBERG_START_MARKER in line:
                found_wiberg = True
                matrix = [[0.0] * n for _ in range(n)]
                continue
                
            if not found_wiberg:
                continue
                
            if line:
                parts = line.split()
                
                if parts[0].lower() == 'atom':
                    current_cols = [int(x) for x in parts[1:]]
                    reading_data = False
                    continue
                
                # Start reading data after separator line
                if all(c in '-' for c in parts[0]):
                    reading_data = True
                    continue
                
                # Process data line
                if reading_data and parts[0].strip().rstrip('.').isdigit():
                    row_idx = int(parts[0].strip().rstrip('.')) - 1
                    values = [float(x) for x in parts[2:]]  # Skip atom type
                    
                    # Store values in matrix
                    for val, col in zip(values, current_cols):
                        matrix[row_idx][col-1] = val
                
            else:  # Empty line
                reading_data = False
                
            if WIBERG_END_MARKER in line:
                break

    if matrix is None:
        raise ValueError("No Wiberg bond order matrix could be generated")
    
    # Convert to numpy array
    matrix_array = np.array(matrix)
    
    # Save to output folder with 4 decimal places formatting
    output_path = output_dir / f'{filename}_wiberg{OUTPUT_EXTENSION}'
    np.savetxt(output_path, matrix_array, delimiter=',', fmt='%.4f')
        
    return matrix_array

def get_wiberg_matrix_size(filename):
    """
    Determines the size of the Wiberg bond index matrix by counting rows.
    """
    input_dir = Path(INPUT_DIR)
    file_path = input_dir / filename
    
    with open(file_path, 'r') as file:
        reading_data = False
        found_wiberg = False
        count = 0
        
        for line in file:
            if WIBERG_START_MARKER in line:
                found_wiberg = True
                continue

            if not found_wiberg:
                continue

            parts = line.strip().split()
            if not parts:  # Empty line
                if reading_data:
                    return count
                continue

            if all(c in '-' for c in parts[0]):  # Separator line
                reading_data = True
                continue

            if reading_data:
                count += 1

    return count