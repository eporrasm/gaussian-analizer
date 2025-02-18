#!/usr/bin/env python3

import argparse
from utils import read_gaussian_matrices, visualize

def main():
    parser = argparse.ArgumentParser(description='Process Gaussian log files to extract matrices.')
    parser.add_argument('filename', help='Name of the log file to process (must be in input folder)')
    
    args = parser.parse_args()
    
    try:
        npa_df, wiberg_df = read_gaussian_matrices(args.filename)
        print(f"Successfully processed {args.filename}")
        print("Files saved in output folder as:")
        print(f"- {args.filename}_npa.csv")
        print(f"- {args.filename}_wiberg.csv")
        visualize(npa_df, wiberg_df)
        print("Graphs created and saved in the output folder.")
    except Exception as e:
        print(f"Error processing file: {e}")

if __name__ == "__main__":
    main()