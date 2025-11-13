# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 22:03:40 2024

@author: tik105
"""

import pandas as pd

# Input and output file paths
input_file = "C://Users//tik105//Desktop//Neonate phenotyping 5-1-2024//Bulk RNA-seq//Null Livers//Results//gene_readcounts_of_all_samples.txt"  # Replace with your text file path
output_file = 'C://Users//tik105//Desktop//Neonate phenotyping 5-1-2024//Bulk RNA-seq//Null Livers//Results//gene_readcounts_of_all_samples.xlsx'  # Replace with your desired Excel file name

# Read the space-delimited text file
try:
    df = pd.read_csv(input_file, delim_whitespace=True)
    # Save to Excel
    df.to_excel(output_file, index=False)
    print(f"File has been saved as {output_file}")
except Exception as e:
    print(f"An error occurred: {e}")
