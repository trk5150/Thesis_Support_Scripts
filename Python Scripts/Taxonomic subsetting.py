# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 15:37:21 2025

@author: tik105
"""

import pandas as pd
import random

# Load the Excel file
file_path = "C://Users//tik105//OneDrive - Harvard University//Atollin Paper//Conservation//BlastPResults_Taxonomy.xlsx"
data = pd.read_excel(file_path)

# Filter rows with "include" set to True
specific_rows = data[data["Include"] == True]

# Function to sample proportional rows from each order
def sample_proportional(data, group_col, fraction=0.2, random_seed=30):
    """
    Randomly samples rows proportionally from each unique group in the specified column.
    :param data: DataFrame containing the data
    :param group_col: Column name to group by
    :param fraction: Fraction of each group to sample
    :param random_seed: Random seed for reproducibility
    :return: DataFrame with sampled rows
    """
    random.seed(random_seed)
    sampled_data = (
        data.groupby(group_col)
        .apply(lambda group: group.sample(frac=fraction, random_state=random_seed))
        .reset_index(drop=True)
    )
    return sampled_data

# Exclude the mandatory rows from the data before sampling
data_for_sampling = data[data["Include"] != True]

# Sample proportional rows from each "order"
fraction_to_sample = 0.1  # Adjust the fraction as needed
sampled_data = sample_proportional(data_for_sampling, group_col="order", fraction=fraction_to_sample)

# Combine sampled rows and specific rows
final_data = pd.concat([sampled_data, specific_rows]).drop_duplicates().reset_index(drop=True)

# Save the resulting data to a new Excel file
output_file = "C://Users//tik105//OneDrive - Harvard University//Atollin Paper//Conservation//BlastPResults_Taxonomy_subsetted.xlsx"
final_data.to_excel(output_file, index=False)

print(f"Proportional sampling completed. The resulting dataset is saved to {output_file}.")
