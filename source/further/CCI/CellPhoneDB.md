# CellPhoneDB

https://www.biorxiv.org/content/10.1101/680926v1.full.pdf

How to use?

First, generate counts matrix (normalized) and metadata (including cell type or functional labeled clusters)

Then, run CellPhoneDB with following commands,

# CellPhoneDB run on linux 
(recommend to use cluster server or workstation)

### We highly recommend using a virtual environment (steps 1 and 2), but can be omitted.

1. Create python > 3.5 virtual-env
$ python -m venv cpdb-venv


2. Activate virtual-env
$ source cpdb-venv/bin/activate


3. Install cellphonedb
$ pip install cellphonedb

# Running with statistical analysis
4. Activate virtual-env if you have not activated previously

$ source cpdb-venv/bin/activate


5. Run in statistical analysis mode using the input files for metadata and counts

$ cellphonedb method statistical_analysis test_meta.txt test_counts.txt

# Result analysis and visualization

### dot plot: indicate activity of ligand-receptor interaction bet. interacting cell types
6. Run after running cellphonedb in either statistical analysis mode or normal mode using the means.csv and pvalues.scv output files

$ cellphonedb plot dot_plot

### Heatmap: indicate frequency of CCIs among various interacting pait-cell types
7. Run after running cellphonedb in either statistical analysis mode or normal mode using the the pvalues.scv output file.

$ cellphonedb plot heatmap_plot meta_data

### Without any optional parameters, you will get total pairs in dot-plot, therefore parsing interacting pairs of geneset is recommended in visualization step. 

8. Validation of CCI

In R, confirm the expression pattern of interacting pairs got from CellPhoneDB process. 
