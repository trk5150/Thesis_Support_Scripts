# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 11:55:41 2024

@author: tik105
"""

import os
import shutil

# Paths
file_list_path = "C://Users//tik105//Desktop//BigScanTest//Prediction vs actual gene lists//Files.txt"  # Replace with the path to your text file
source_directory = "C://Users//tik105//Desktop//LINCs sets//GMT parsed only small molecules"  # Replace with your source directory
target_directory = "C://Users//tik105//Desktop//BigScanTest//Prediction vs actual gene lists"  # Replace with your target directory

# Ensure the target directory exists
os.makedirs(target_directory, exist_ok=True)

# Read the file names from the text file
with open(file_list_path, 'r') as file:
    file_names = [line.strip() for line in file if line.strip()]

# Search and copy files
for file_name in file_names:
    source_file_path = os.path.join(source_directory, file_name)
    if os.path.isfile(source_file_path):
        shutil.copy(source_file_path, target_directory)
        print(f"Copied: {file_name}")
    else:
        print(f"File not found: {file_name}")

# Paths to the files containing gene names
file1_path = "C://Users//tik105//Desktop//BigScanTest//Prediction vs actual gene lists//Epirizole.txt"  # Replace with the path to the first file
#file1_path = "C://Users//tik105//Desktop//BigScanTest//Prediction vs actual gene lists//Moxifloxacin.txt"  # Replace with the path to the first file

file2_path = "C://Users//tik105//Desktop//BigScanTest//Prediction vs actual gene lists//LCP001_MCF10A.WT_24H_L20_epirizole_2.5uM up.txt"  # Replace with the path to the second file
#file2_path = "C://Users//tik105//Desktop//BigScanTest//Prediction vs actual gene lists//ASG002_JHH5_24H_L13_moxifloxacin_10uM up.txt"
#file2_path = "C://Users//tik105//Desktop//BigScanTest//Prediction vs actual gene lists//REP.B018_A375_24H_N17_moxifloxacin_0.03uM up.txt"

# Read the gene names from both files, ignoring lines with only "-"
def read_gene_names(file_path):
    with open(file_path, 'r') as file:
        return {line.strip() for line in file if line.strip() and line.strip() != "-"}

# Load gene names from both files
genes_file1 = read_gene_names(file1_path)
genes_file2 = read_gene_names(file2_path)

# Find the common gene names
common_genes = genes_file1 & genes_file2  # Intersection of the sets

# Output the result
print(f"Number of common genes: {len(common_genes)}")