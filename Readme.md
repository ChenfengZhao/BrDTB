# BRDTB
## Overview
Add a few sentence to discribe this project.

## Requirements
- Recommended OS: Linux (e.g., Ubuntu)
- Python3 (e.g., 3.12.3)
- Pip3
- requirements.txt

## Installation Guide

Run the following command to install the dependencies in a virtual environment
```
python3 -m venv .venv # create a virtual environment
source .venv/bin/activate # start the virtual environment
pip install -r requirements.txt # install python dependencies
```

Run the following command to exit the virtual environment once the program is finished
```
deactivate
```

## Usage Example
1. Download the source code of brdtb from Github.

2. Prepare excel and mzXML files in ./data/ (e.g., test_peptide_BrDTB.xlsx, 20240802_BBP_Peptide_1_BrDTB.mzXML).

3. Run the program using the following command
```
python3 brdtb_cal.py
```

4. Results will be generated in ./results/. If users want to rerun the same peptide, delete the following row in ./results/result_summary.xlsx. Otherwise, the program will not process the peptides that have already been processed.


## License
[Apache_2.0_license]: http://www.apache.org/licenses/LICENSE-2.0

The source code of this project is released under the [Apache 2.0 License][Apache_2.0_license].

## Citation
If you think BRDTB is helpful for your research, please cite the following paper: