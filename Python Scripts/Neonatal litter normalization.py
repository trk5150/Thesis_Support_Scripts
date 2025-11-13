# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 16:22:02 2025

@author: tik105
"""

import pandas as pd

def process_genotype_data(file_path):
    # Load the Excel file
    xls = pd.ExcelFile(file_path)
    
    # Read the sheet named "genotype sorted"
    df = xls.parse("Genotype sorted")
    
    # Ensure correct column names and types
    numeric_cols = ["Mass (g)", "Blood glucose (mg/dL)", "Ketones (mg/dL)", "Lactate"]
    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors='coerce')
    
    # Calculate average values for each Litter
    litter_avg = df.groupby("Litter")[numeric_cols].mean().rename(columns=lambda x: f"Avg {x}")
    
    # Merge the averages back into the main dataframe
    df = df.merge(litter_avg, on="Litter", how="left")
    
    # Calculate normalized values as percent difference from litter average
    for col in numeric_cols:
        df[f"Normalized {col}"] = ((df[col] - df[f"Avg {col}"]) / df[f"Avg {col}"]) * 100
    
    # Save the modified dataframe to a new Excel file
    output_file = "C://Users//tik105//Desktop//Neonate phenotyping 5-1-2024//processed_genotype_data.xlsx"
    df.to_excel(output_file, index=False)
    
    print(f"Processed data saved to {output_file}")
    return df

# Example usage
file_path = "C://Users//tik105//Desktop//Neonate phenotyping 5-1-2024//Mass, BG and ketones.xlsx"  # Replace with your actual file path
df_processed = process_genotype_data(file_path)
