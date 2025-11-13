# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 13:48:32 2024

@author: tik105
"""

import pandas as pd

# Define the file path and threshold values
file_path = "C://Users//tik105//Desktop//BigScanTest//Map Consistency Output//5050 full sc-b//SOM_Analysis.txt"  # Update this to the path of your tab-delimited file
structure_threshold = .0750  # Set your desired threshold for 'Structure'
quality_threshold = .057  # Set your desired threshold for 'quality'

# Read the tab-delimited file into a pandas DataFrame
df = pd.read_csv(file_path, sep='\t')

# Count rows where 'Structure' exceeds the threshold
structure_count = len(df[df['Structure'] > structure_threshold])

# Count rows where 'quality' exceeds the threshold
quality_count = len(df[df['quality'] > quality_threshold])

count = len(df[(df['Structure'] > structure_threshold) & (df['quality'] > quality_threshold)])

# Print the counts
print(f"Number of rows with 'Structure' > {structure_threshold}: {structure_count}")
print(f"Number of rows with 'quality' > {quality_threshold}: {quality_count}")
print(f"Number of rows with 'Structure' > {structure_threshold} and 'quality' > {quality_threshold}: {count}")
