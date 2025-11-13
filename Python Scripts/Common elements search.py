# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 14:11:06 2024

@author: tik105
"""

import os
import pandas as pd

# Define the paths


target_file_path = "C://Users//tik105//Desktop//BigScanTest//Prediction vs actual gene lists//Moxifloxacin.txt"
directory_path = "C://Users//tik105//Desktop//LINCs sets//GMT parsed only small molecules"  # Replace with your source directory
#file1_path = "C://Users//tik105//Desktop//BigScanTest//Prediction vs actual gene lists//Epirizole.txt"  # Replace with the path to the first file
#file1_path = "C://Users//tik105//Desktop//BigScanTest//Prediction vs actual gene lists//Moxifloxacin.txt"  # Replace with the path to the first file


# Load the target gene list
with open(target_file_path, 'r') as file:
    target_genes = set(file.read().splitlines())

# Prepare a list to store results
results = []

# Iterate through each file in the directory
for file_name in os.listdir(directory_path):
    file_path = os.path.join(directory_path, file_name)
    if os.path.isfile(file_path):  # Check if it's a file
        # Load the genes from the current file
        with open(file_path, 'r') as file:
            genes = set(file.read().splitlines())
        
        # Find the common genes
        common_genes = target_genes.intersection(genes)
        
        # Append the results
        results.append({'File Name': file_name, 'Common Genes': len(common_genes)})

# Create a DataFrame to store the results
results_df = pd.DataFrame(results)

# Calculate statistics
statistics = {
    "Average": results_df["Common Genes"].mean(),
    "Standard Deviation": results_df["Common Genes"].std(),
    "Minimum": results_df["Common Genes"].min(),
    "Maximum": results_df["Common Genes"].max()
}

# Convert statistics to a DataFrame for better organization
statistics_df = pd.DataFrame(statistics, index=["Value"])

print("\nStatistics:")
print(statistics_df)

# Display the DataFrame
#print(results_df)

# Save the results to a CSV if needed
results_df.to_csv("C://Users//tik105//Desktop//BigScanTest//Prediction vs actual gene lists//Moxifloxacin_common_genes_summary.csv", index=False)
statistics_df.to_csv("C://Users//tik105//Desktop//BigScanTest//Prediction vs actual gene lists//Moxifloxacin_common_genes_statistics.csv", index=False)
