# CellPhoneDB

https://www.biorxiv.org/content/10.1101/680926v1.full.pdf

How to use?

First, generate counts matrix (normalized) and metadata (including cell type or functional labeled clusters)

Then, run CellPhoneDB with following commands,

# CellPhoneDB run on linux (recommend to use cluster server or workstation)
We highly recommend using a virtual environment (steps 1 and 2), but can be omitted.

1. Create python > 3.5 virtual-env
python -m venv cpdb-venv


2. Activate virtual-env
source cpdb-venv/bin/activate


3. Install cellphonedb
pip install cellphonedb

# Running with statistical analysis
4. Activate virtual-env if you have not activated previously
source cpdb-venv/bin/activate


5. Run in statistical analysis mode using the input files for metadata and counts
cellphonedb method statistical_analysis test_meta.txt test_counts.txt
