# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 17:21:58 2025

@author: tik105
"""

import pandas as pd

# Load CSV and Excel files
csv_file = "C://Users//tik105//Downloads//ENSMUSG00000070577_correlation.csv"  # Replace with actual file path
excel_file = "C://Users//tik105//Desktop//Neonate phenotyping 5-1-2024//Bulk RNA-seq//Null pancreas//Results//DiffEx.xlsx"  # Replace with actual file path

# Read the CSV and Excel files
df_csv = pd.read_csv(csv_file)
df_excel = pd.read_excel(excel_file)

# Merge the dataframes on 'symbol'
df_merged = df_csv.merge(df_excel, on="symbol", how="inner")

# Filter where padj > 0.05
df_filtered = df_merged[df_merged["padj"] < 0.05]

# Add 'Agree' column based on conditions
df_filtered["Agree"] = ((df_filtered["correlation"] > 0) & (df_filtered["log2FoldChange"] < 0)) | \
                        ((df_filtered["correlation"] < 0) & (df_filtered["log2FoldChange"] > 0))
df_filtered["Agree"] = df_filtered["Agree"].astype(int)

# Select relevant columns
result_df = df_filtered[["symbol", "correlation", "padj", "log2FoldChange", "Agree"]]

# Save the result to a new Excel file
result_df.to_excel("C://Users//tik105//Desktop//Neonate phenotyping 5-1-2024//Bulk RNA-seq//Null pancreas//Results//attie_filtered_data.xlsx", index=False)
