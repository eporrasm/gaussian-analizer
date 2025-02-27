#!/usr/bin/env python3

from pathlib import Path
from utils import read_gaussian_matrices, visualize
from config import INPUT_DIR

def get_input_file(prompt):
    while True:
        filename = input(prompt)
        if Path(INPUT_DIR).joinpath(filename).is_file():
            return filename
        print(f"Error: File not found in {INPUT_DIR} directory. Please try again.")

def main():
    print("Please provide the following Gaussian log files from the input folder:")
    
    # Get the three input files
    first_fragment = get_input_file("Enter filename for first fragment: ")
    other_fragments = get_input_file("Enter filename for repeating fragments: ")
    final_state = get_input_file("Enter filename for final state after simulation: ")
    
    try:
        # Process first fragment
        npa_first, wiberg_first = read_gaussian_matrices(first_fragment)
        print(f"\nSuccessfully processed {first_fragment}")
        print(f"Files saved as {first_fragment}_npa.csv and {first_fragment}_wiberg.csv")
        
        # Process other fragments
        npa_other, wiberg_other = read_gaussian_matrices(other_fragments)
        print(f"\nSuccessfully processed {other_fragments}")
        print(f"Files saved as {other_fragments}_npa.csv and {other_fragments}_wiberg.csv")
        
        # Process final state and create visualizations
        npa_final, wiberg_final = read_gaussian_matrices(final_state)
        print(f"\nSuccessfully processed {final_state}")
        print(f"Files saved as {final_state}_npa.csv and {final_state}_wiberg.csv")
        
        # Create visualizations for the final state
        visualize(npa_final, wiberg_final)
        print("\nGraphs created and saved in the output folder.")
        
    except Exception as e:
        print(f"Error processing files: {e}")

if __name__ == "__main__":
    main()