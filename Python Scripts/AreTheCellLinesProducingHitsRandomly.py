# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 17:24:53 2024

@author: tik105
"""

from scipy.stats import hypergeom
import pandas as pd

# Load the Excel file into a DataFrame
excel_file_path = "C://Users//tik105//OneDrive - Harvard University//SOMSCAN draft//SOMSCAN figures//Figures//Figure 4 l1000 screening//CellLineSigsandHits.xlsx"  # Replace with the actual path to your file
data = pd.read_excel(excel_file_path)

# Create empty list to store results
overrepresented_cell_lines = []

# Total number of tests and hits in the dataset
total_tests = data['Signatures Tested'].sum()
total_hits = data['Hits'].sum()

# Loop through each cell line and perform the hypergeometric test
for index, row in data.iterrows():
    cell_line = row['Cell Line']
    cell_line_tests = row['Signatures Tested']
    cell_line_hits = row['Hits']
    
    # Perform the hypergeometric test
    # N: total_tests, K: total_hits, n: cell_line_tests, k: cell_line_hits
    p_value = hypergeom.sf(cell_line_hits - 1, total_tests, total_hits, cell_line_tests)  # sf is for P(X >= k)
    
    # If the p-value is significant, add the cell line to the results
    if p_value < 0.05:  # Adjust this threshold based on your desired significance level
        overrepresented_cell_lines.append((cell_line, p_value))

# Print the results
print("Overrepresented Cell Lines (p < 0.05):")
for cell_line, p_value in overrepresented_cell_lines:
    print(f"{cell_line}: p-value = {p_value:.5f}")