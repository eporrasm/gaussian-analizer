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
        print('The fragment 1 is related to the carbonized structure')
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

def visualize(npa_df, wiberg_df, npa_first, wiberg_first, npa_other, wiberg_other):
    """Creates visualizations of the NPA and Wiberg bond index matrices."""
    pd.options.mode.chained_assignment = None
    fragment_division = user_input(len(npa_df.index))
    graph_wiberg(wiberg_df, fragment_division)
    print("Wiberg heatmap completed. Files created.\n")
    graph_npa(npa_df, npa_first, npa_other,  fragment_division)
    print("NPA completed. Files created.\n")
   


def graph_npa(df, npa_first, npa_other, fragment_division):
    """Creates a horizontal bar plot of the Natural Population Analysis (NPA) data.
    Parameters:
    df (pd.DataFrame): The DataFrame containing the NPA data.
    npa_first: The DataFrame containing the NPA data of the first fragment (carbonized).
    npa_other: The DataFrame containing the NPA data of the other fragments.
    fragment_division (list): A list of integers representing the division of the fragments.
    Returns:
    None: The function saves the plot as a PNG file and the data as a CSV file.
    """
    row = df
    counter = 1
    row_split = list()

    df['Atom_Num'] = df['Atom_Num'].astype('int') 
    npa_first['Atom_Num'] = npa_first['Atom_Num'].astype('int')
    npa_other['Atom_Num'] = npa_other['Atom_Num'].astype('int')
    npa_first['Natural_Charge_Before'] = npa_first['Natural_Charge'].astype('float')
    npa_other['Natural_Charge_Before'] = npa_other['Natural_Charge'].astype('float')
    npa_first = npa_first[['Atom_Num', 'Natural_Charge_Before']]
    npa_other = npa_other[['Atom_Num', 'Natural_Charge_Before']]         

    for i in range(len(fragment_division) - 1):
        if i == 0:
            npa_data = npa_first
        else:
            npa_data = npa_other
            npa_data['Atom_Num'] += (fragment_division[i]) - (fragment_division[i-1])

        new_row = row.loc[fragment_division[i]:fragment_division[i+1]-1]
        new_row = new_row.join(npa_data.set_index('Atom_Num'), on='Atom_Num')
        new_row['Natural_Charge_Diff'] = new_row['Natural_Charge_Before'] - new_row['Natural_Charge']
        row_split.append(new_row)

    filter_type, filter_value = input_filter_npa()
    for i in range(len(row_split)):
        row_split[i] = filter_npa(row_split[i], filter_type, filter_value)

    orderby = order_by_npa()
    
    for row_data in row_split:
        if row_data.empty:
            print(f'Warning: The fragment {counter} is empty! No image will be generated.')
            counter += 1
            continue
    
        if orderby == 'charge':
            row_data['Natural_Charge_Diff_Abs'] = row_data['Natural_Charge_Diff'].abs()
            row_data = row_data.sort_values(by='Natural_Charge_Diff_Abs', ascending=True)

        # Determine the figure size based on the data size
        if row_data.shape[0] > 50 or row_data.shape[1] > 60:
            plt.figure(figsize=(10, 40))
        elif row_data.shape[0] > 30 or row_data.shape[1] > 30:
            plt.figure(figsize=(10, 30))
        else:
            plt.figure(figsize=(10, 20))

        # Create an array for the y positions
        y_positions = np.arange(len(row_data['Atom_Num']))

        # Define the width of the bars
        bar_width = 0.35

        # Create the horizontal bar plot with slight offset
        bars = plt.barh(y_positions - bar_width/2, row_data['Natural_Charge'], height=bar_width, align='center', color='skyblue', label='After')
        bars2 = plt.barh(y_positions + bar_width/2, row_data['Natural_Charge_Before'], height=bar_width, align='center', color='yellow', label='Before')

        # Add the name of each entry at the end of the bar
        for bar, atom_num, y in zip(bars, row_data['Atom_Num'], y_positions):
            x = bar.get_width()  # Width of the bar (value)
            plt.text(x, y - bar_width/2, f'Atom {atom_num} (after): {x:.2f}', ha='left', va='center', color='black')

        for bar, atom_num, y in zip(bars2, row_data['Atom_Num'], y_positions):
            x = bar.get_width()  # Width of the bar (value)
            plt.text(x, y + bar_width/2, f'Atom {atom_num} (before): {x:.2f}', ha='left', va='center', color='black')

        # Add labels and title
        plt.ylabel('Atoms (identified by number)')
        plt.xlabel('Natural Charge')
        plt.title('Horizontal Bar Diagram of Natural Charges')
        plt.yticks(y_positions, row_data['Atom_Num'])  # Set y-ticks to Atom_Num

        # Add extras
        plt.grid(visible=True, axis='x')  # Add grid lines only for the x-axis
        plt.legend()  # Add legend to distinguish 'Before' and 'After'

        # Show the plot
        plt.tight_layout()  # Adjust layout to prevent clipping of labels
        plt.savefig(f'./output/npa/npa_fragment{counter}.png', dpi=300, bbox_inches='tight')
        # Save the plot's data as CSV
        row_data[['Atom_Num', 'Natural_Charge_Before', 'Natural_Charge', 'Natural_Charge_Diff']].to_csv(
            f'./output/npa/logs/npa_fragment{counter}.csv',
            index=False)
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
    """Filter the NPA data based on the absolute value of the charge change."""
    return df[df['Natural_Charge_Diff'].abs() > filter_value]

def order_by_npa():
    """Ask the user if they want to order the NPA data by the atom number or the charge value.
    Returns:
    str: The column name to use for ordering."""
    order_by = input('Order the NPA data by atom number or charge value change? ((1)ATOM / (2)charge) ')
    if order_by.lower() in ['atom', 'charge', '', '1', '2']:
        if order_by.lower() in ['atom', '1', '']:
            return 'atom'
        else:
            return 'charge'
    else:
        print('Invalid answer\n')
        return order_by_npa()

def graph_wiberg(df, fragment_division):
    """Creates a heatmap of the Wiberg bond index matrix.
    Parameters:
    df (np.ndarray): The Wiberg bond index matrix.
    fragment_division (list): A list of integers representing the division of the fragments.
    returns:
    None: The function saves the plot as a PNG file and the data as a CSV file.
    """
    matrixes = list()
    atom_list = list()

    # Split the matrix into its fragments
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

        # Save the plot
        plt.savefig(f'./output/wiberg/wiberg_fragment{i+2}xfragment1.png', dpi=300, bbox_inches='tight')
        # Save the plot's data as CSV
        np.savetxt(f'./output/wiberg/logs/wiberg_fragment{i+2}xfragment1.csv', matrix, delimiter=',', fmt='%.4f')
        # Save the atom labels as CSV
        pd.DataFrame({'Fragment 1': x_labels}).to_csv(
            f'./output/wiberg/logs/wiberg_fragment{i+2}xfragment1_xlabels.csv', index=False)
        pd.DataFrame({f'Fragment {i+2}': y_labels}).to_csv(
            f'./output/wiberg/logs/wiberg_fragment{i+2}xfragment1_ylabels.csv', index=False)

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
        atoms1: The list of atoms for the first fragment.
        atoms2: The list of atoms of the fragment compared to the first fragment.
        threshold (float): The threshold value.
    
    Returns:
        numpy.ndarray: The modified matrix.
        list: The modified list of atoms for the first fragment.
        list: The modified list of atoms for the other fragment. 
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
