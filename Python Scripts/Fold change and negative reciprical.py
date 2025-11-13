# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 14:00:01 2024

@author: tik105
"""

import pandas as pd
import numpy as np

# Load RNA-seq count data from Excel
file_path = "D://Bulk RNA-seq results//24231-01-Analysis-09262024_142004//DE calcs.xlsx"  # Replace with your Excel file path

# Load the data
df = pd.read_excel(file_path)
# Controls are in "DMSO 1" and "DMSO 2"
control_cols = ['DMSO 1', 'DMSO 2']

# List of treatment column pairs, e.g., ['Treatment_1_1', 'Treatment_1_2', 'Treatment_2_1', 'Treatment_2_2']
# Add your treatment pairs here, for example:
treatment_cols = [['Colch 1', 'Colch 2'], ['Epirizole 1', 'Epirizole 2'], ["LY-294002 1", "LY-294002 2"], ["Moxifloxacin 1", "Moxifloxacin 2"], ["Rottlerin 1", "Rottlerin 2"], ["Vardenafil 1", "Vardenafil 2"], ["Rimon 1", "Rimon 2"]]  # Add more if needed

# Calculate the mean of the control group
df['Mean_Control'] = df[control_cols].mean(axis=1)

# Loop through each treatment pair and calculate fold change
for treatment_pair in treatment_cols:
    # Calculate the mean of the treatment group
    mean_treatment_col = f'Mean_{treatment_pair[0].split("_")[0]}'  # e.g., "Mean_Treatment_1"
    df[mean_treatment_col] = df[treatment_pair].mean(axis=1)
    
    # Calculate fold change
    fold_change_col = f'Fold_Change_{treatment_pair[0].split("_")[0]}'
    df[fold_change_col] = df[mean_treatment_col] / df['Mean_Control']
    
    # Apply negative reciprocal for values less than 1
    #df[fold_change_col] = df[fold_change_col].apply(lambda x: -1/x if x < 1 else x)
    
    # Calculate log2 fold change
    log2_fold_change_col = f'Log2_Fold_Change_{treatment_pair[0].split("_")[0]}'
    df[log2_fold_change_col] = np.log2(df[fold_change_col])
    
# Create two new rows: one for values > 0 and one for values < 0
row_greater_than_zero = {col: (df[col] > 1).sum() for col in df.columns if col.startswith('Fold_Change')}
row_less_than_zero = {col: (df[col] < 1).sum() for col in df.columns if col.startswith('Fold_Change')}

# Add 'Gene' label to the new rows
row_greater_than_zero['Gene'] = 'Values > 0'
row_less_than_zero['Gene'] = 'Values < 0'

# Append these new rows to the DataFrame
df = df.append([row_greater_than_zero, row_less_than_zero], ignore_index=True)

# Optionally save the output to a new Excel file
output_file = "D://Bulk RNA-seq results//24231-01-Analysis-09262024_142004//DE fold calcs.xlsx"
df.to_excel(output_file, index=False)
print(f"Results saved to {output_file}")