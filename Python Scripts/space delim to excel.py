# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 16:47:23 2024

@author: tik105
"""

import pandas as pd

# Specify the input and output file paths
input_file =   "C://Users//tik105//Desktop//Neonate phenotyping 5-1-2024//Bulk RNA-seq//Null pancreas//gene_readcounts_of_all_samples.txt"# replace with your space-delimited file
output_file = "C://Users//tik105//Desktop//Neonate phenotyping 5-1-2024//Bulk RNA-seq//Null pancreas//Counts.xlsx"

# Read the space-delimited file into a DataFrame
df = pd.read_csv(input_file, delim_whitespace=True)

# Export the DataFrame to an Excel file
df.to_excel(output_file, index=False)

print(f"Data has been successfully exported to {output_file}")
