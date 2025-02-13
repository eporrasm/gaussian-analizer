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

1. Place your Gaussian .log files in the `input` folder
2. Run the analyzer:
```bash
python gaussian-analizer.py your_file.log
```

The program will generate two CSV files in the `output` folder:
- `your_file.log_npa.csv`: Natural Population Analysis matrix
- `your_file.log_wiberg.csv`: Wiberg bond index matrix

## Example
```bash
# Copy your Gaussian log file to input folder
cp path/to/your/file.log input/

# Run the analyzer
python gaussian-analizer.py file.log
```
