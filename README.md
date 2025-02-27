# Gaussian Log File Analyzer

A command-line tool to extract and analyze Natural Population Analysis and Wiberg bond index matrices from Gaussian log files.

## Requirements
- Python 3.7 or higher

## Installation

1. Clone this repository

2. Create and activate a virtual environment:
```bash
# Create virtual environment
python -m venv venv

# Activate it on Linux/Mac
source venv/bin/activate

# Or on Windows
# venv\Scripts\activate
```

3. Install the required dependencies:
```bash
# Make sure you're in the project directory and your virtual environment is activated
pip install -r requirements.txt
```

## Usage

1. Place your Gaussian .log files in the `input` folder. You need three files:
   - First fragment log file
   - Repeating fragments log file
   - Final state log file (after simulation)

2. Run the analyzer:
```bash
python gaussian-analizer.py
```

3. When prompted, enter the filenames for:
   - First fragment
   - Repeating fragments
   - Final state

The program will generate CSV files in the `output` folder for each input file:
- `*_npa.csv`: Natural Population Analysis matrix
- `*_wiberg.csv`: Wiberg bond index matrix

Then the program will prompt for fragment information and generate visualizations for the final state data.

## Example
```bash
# Copy your Gaussian log files to input folder
cp path/to/first_fragment.log input/
cp path/to/other_fragments.log input/
cp path/to/final_state.log input/

# Run the analyzer and follow the prompts
python gaussian-analizer.py
```
