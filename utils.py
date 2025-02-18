import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from config import INPUT_DIR, OUTPUT_DIR, OUTPUT_EXTENSION, NPA_START_MARKER, NPA_END_MARKER, WIBERG_START_MARKER, \
        WIBERG_END_MARKER, NPA_COLUMNS, YES, NO

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

def user_input(max_atom_num):
    """Ask the user for the number of fragments and the number of atoms in each fragment.
    parameters:
    max_atom_num (int): The total number of atoms in the molecule.
    Returns:
    list: A list containing the division of the fragments."""
    print(f'The model has a total of {max_atom_num} atoms')
    max_atom_num = int(max_atom_num)

    while True:
        fragments = list()
        division = [0]
        atom_num = 0

        print('Number of fragments?')
        n_fragments = int(input())

        for i in range(n_fragments-1):
            fragments.append(int(input(f"How many atoms does the fragment {i+1} have? ")))

        if sum(fragments) >= max_atom_num:
            raise ValueError("The sum of the fragments is greater than the total number of atoms")
        
        print(f"Then the last fragment will have {max_atom_num - sum(fragments)} atoms")
        fragments.append(max_atom_num - sum(fragments))
        str_fragments = ""

        for i in range(n_fragments):
            str_fragments += f"The fragment {i+1} goes from {int(atom_num) + 1} "
            atom_num += fragments[i]
            str_fragments += f"to {atom_num}\n"   
            division.append(atom_num)

        while True:
            str_fragments += f"Is this correct? Y/n "
            correct = input(str_fragments)
            if correct in YES:
                return division
            elif correct in NO:
                break
            else:
                print('Invalid answer\n')

def visualize(npa_df, wiberg_df):
    """Creates visualizations of the NPA and Wiberg bond index matrices."""
    fragment_division = user_input(len(npa_df.index))
    graph_wiberg(wiberg_df, fragment_division)
    print("Wiberg heatmap created and saved in the output folder.")
    graph_npa(npa_df, fragment_division)
   


def graph_npa(df, fragment_division):
    """Creates a horizontal bar plot of the Natural Population Analysis (NPA) data.
    Parameters:
    df (pd.DataFrame): The DataFrame containing the NPA data.
    fragment_division (list): A list of integers representing the division of the fragments.
    """
    row = df
    counter = 1
    row_split = list()

    for i in range(len(fragment_division) - 2):
        row_split.append(row.loc[fragment_division[i]:fragment_division[i+1]-1])

    i += 1
    row_split.append(row.loc[fragment_division[i]:fragment_division[i+1]])

    filter_type, filter_value = input_filter_npa()
    for i in range(len(row_split)):
        row_split[i] = filter_npa(row_split[i], filter_type, filter_value)

    for row_data in row_split:
        if row_data.empty:
            print(f'Warning: The fragment {counter} is empty! No image will be generated.')
            counter += 1
            continue
        # Create the horizontal bar plot
        plt.figure(figsize=(10, 30))  # Adjust the figure size as needed
        bars = plt.barh(row_data['Atom_Num'], row_data['Natural_Charge'], align='center', color='skyblue')

        # Add the name of each entry at the end of the bar
        for bar, atom_num in zip(bars, row_data['Atom_Num']):
            # Get the x and y positions for the label
            x = bar.get_width()  # Width of the bar (value)
            y = bar.get_y() + bar.get_height() / 2  # Center of the bar vertically
            
            # Add the label in the format "Atom {Atom_Num}: {Natural_Charge}"
            plt.text(x, y, f'Atom {atom_num}: {x:.2f}', ha='left', va='center', color='black')

        # Add labels and title
        plt.ylabel('Categories')
        plt.xlabel('Natural Charge')
        plt.title('Horizontal Bar Diagram of Natural Charges')
        plt.yticks([])  # Hide default y-axis labels since we're adding custom ones

        # Add extras
        plt.grid(visible=True, axis='x')  # Add grid lines only for the x-axis

        # Show the plot
        plt.tight_layout()  # Adjust layout to prevent clipping of labels
        plt.savefig(f'./output/npa_fragment{counter}.png', dpi=300, bbox_inches='tight')
        counter += 1

def input_filter_npa():
    """Ask the user if they want to filter the NPA data based on the charge value.
    Returns:
    tuple: A tuple containing the filter type and the filter value."""
    filter = input('Only show atoms of NPA with a charge greater than a certain value? Y/n ')
    if filter in YES:
        filter_type = 'abs_value'
        filter_value = input('Enter the value (default 0.3): ')
        if filter_value == '':
            filter_value = 0.3
        return filter_type, float(filter_value)
    elif filter in NO:
        return 'abs_value', 0
    else:
        print('Invalid answer\n')
        return input_filter_npa()

def filter_npa(df, filter_type, filter_value):
    """Filter the NPA data based on the specified filter type and value.
    Parameters:
    df (pd.DataFrame): The DataFrame containing the NPA data.
    filter_type (str): The type of filter to apply.
    filter_value (float): The value to use for the filter.
    Returns:
    pd.DataFrame: The filtered DataFrame."""

    match filter_type:
        case 'abs_value':
            return abs_value_npa(df, filter_value)
        case _:
            pass

def abs_value_npa(df, filter_value):
    """Filter the NPA data based on the absolute value of the charge."""
    return df[df['Natural_Charge'].abs() > filter_value]

def graph_wiberg(df, fragment_division):
    """Creates a heatmap of the Wiberg bond index matrix.
    Parameters:
    df (np.ndarray): The Wiberg bond index matrix.
    fragment_division (list): A list of integers representing the division of the fragments.
    """
    matrixes = list()
    atom_list = list()

    for i in range(1, len(fragment_division)-1):
        matrix = list()
        for j in range(fragment_division[1]):
            matrix.append(df[j][fragment_division[i]:fragment_division[i+1]])
        matrixes.append(matrix)


    atoms_1 = list(range(1, fragment_division[1] + 1))
    tol = input_filter_wiberg()
    for i in range(len(matrixes)):
        atoms_2 = list(range(fragment_division[i+1] + 1, fragment_division[i+2] + 1))
        matrixes[i], a1, a2 = edit_matrix(matrixes[i], atoms_1, atoms_2, tol)
        atom_list.append((a2, a1))
    
    for i, matrix in enumerate(matrixes):
        # Extract the labels for the current matrix
        x_labels = atom_list[i][0]  # Labels for the x-axis
        y_labels = atom_list[i][1]  # Labels for the y-axis
        if len(y_labels) > 100 or len(x_labels) > 100:
            figsize = (80, 80)
        elif len(y_labels) > 50 or len(x_labels) > 50:
            figsize = (50, 50)
        else:
            figsize = (15, 15)
        if matrix.size == 0:
            print(f'Warning: The matrix fragment 1 x fragment {i+2} is empty! No image will be generated.')
            continue
        # Create the heatmap
        plt.figure(figsize=figsize) 
        ax = sns.heatmap(matrix, annot=True, cmap='viridis',fmt=".2f", xticklabels=x_labels, yticklabels=y_labels)

        # Add labels and title
        plt.xlabel(f'Fragment {i+2} (Atom #)')
        plt.ylabel('Fragment 1 (Atom #)')
        plt.title(f'Heatmap of Fragment {i+2} x Fragment 1')
        
        # Rotate x-axis labels for better readability
        plt.xticks(rotation=45)
        plt.yticks(rotation=0)

        # Show the plot
        plt.savefig(f'./output/wiberg_fragment{i+2}xfragment1.png', dpi=300, bbox_inches='tight')

def input_filter_wiberg():
    """Ask the user if they want to filter the Wiberg bond index matrix based on the value.
    Returns:
    float: The threshold value to use for filtering."""
    filter = input('Only show atoms with a Wiberg bond index greater than a certain value? Y/n ')
    if filter in YES:
        filter_value = input('Enter the value (default 0.1): ')
        if filter_value == '':
            filter_value = 0.1
        return float(filter_value)
    elif filter in NO:
        return 0
    else:
        print('Invalid answer\n')
        return input_filter_npa()
    
def edit_matrix(matrix, atoms_1, atoms_2, threshold):
    """
    Edit the matrix by deleting columns and rows where all values are below the threshold.
    
    Parameters:
        matrix: The input matrix.
        threshold (float): The threshold value.
    
    Returns:
        numpy.ndarray: The modified matrix.
    """
    matrix = np.array(matrix)
    # Delete columns where all values are below the threshold
    col_mask = np.any(matrix >= threshold, axis=0)  # Keep columns where at least one value >= threshold
    matrix = matrix[:, col_mask]

    # We keep the values of atoms that are true in col_mask
    atoms_2 = [atom for atom, keep in zip(atoms_2, col_mask) if keep]
    
    # Delete rows where all values are below the threshold
    row_mask = np.any(matrix >= threshold, axis=1)  # Keep rows where at least one value >= threshold
    matrix = matrix[row_mask, :]
    atoms_1 = [atom for atom, keep in zip(atoms_1, row_mask) if keep]
    
    
    return matrix, atoms_1, atoms_2
