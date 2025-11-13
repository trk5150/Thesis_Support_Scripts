# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 17:34:11 2024

@author: tik105
"""
import math
import numpy as np
import pandas as pd
from scipy.stats import kruskal
import scikit_posthocs as sp

# Replace 'file_path.xlsx' with your actual file path
file_path = 'C://Users//tik105//Desktop//GTT false positive estimation//3-23-2018.xlsx'
df = pd.read_excel(file_path)





def split_dataframe_randomly(df, n_groups):
    """
    Splits a DataFrame into a set number of random groups.
    
    Parameters:
    df (pd.DataFrame): The DataFrame to be split.
    n_groups (int): The number of groups to split the DataFrame into.
    
    Returns:
    list of pd.DataFrame: A list containing the DataFrames for each group.
    """
    # Number of rows in the DataFrame
    n_rows = len(df)
    
    # Shuffle the DataFrame indices
    shuffled_indices = np.random.permutation(n_rows)
    
    # Determine the size of each group
    group_sizes = np.array_split(shuffled_indices, n_groups)
    
    # Split the DataFrame into groups
    groups = [df.iloc[indices].reset_index(drop=True) for indices in group_sizes]
    
    return groups

def calculate_auc_per_row(df):
    """
    Calculates the area under the curve (AUC) for each row in the DataFrame using the trapezoidal rule.
    The time points are extracted from the column names.
    
    Parameters:
    df (pd.DataFrame): The input DataFrame containing treatment and time course data.
    
    Returns:
    pd.Series: A series containing the AUC values for each row.
    """
    # Extract time points from column names (all except the first column which is 'Treatment')
    time_points = df.columns[1:].astype(float).values
    
    # Extract the time course data (all except the first column, which is 'Treatment')
    time_course_data = df.iloc[:, 1:].values
    
    # Calculate AUC for each row using the trapezoidal rule
    auc_values = [np.trapz(y=row, x=time_points) for row in time_course_data]
    
    # Return the AUC values as a pandas Series
    return pd.Series(auc_values, index=df.index)

def replace_labels_and_combine(dfs, base_label='Treatment'):
    """
    Replaces labels in the first column of each DataFrame with a unique identifier
    and combines all DataFrames into a single DataFrame.
    
    Parameters:
    dfs (list of pd.DataFrame): List of DataFrames to be processed and combined.
    base_label (str): Base label to use for unique identifiers (default is 'Group').
    
    Returns:
    pd.DataFrame: Combined DataFrame with updated labels.
    """
    combined_df = pd.DataFrame()
    
    for i, df in enumerate(dfs):
        # Create a unique identifier for the current DataFrame
        unique_identifier = f"{base_label}{i+1}"
        
        # Replace labels in the first column with the unique identifier
        df[df.columns[0]] = unique_identifier
        
        # Append to the combined DataFrame
        combined_df = pd.concat([combined_df, df], ignore_index=True)
    
    return combined_df

def perform_dunns_test(df, group_column, auc_column):
    """
    Performs Dunn's multiple comparison test on AUC values after a Kruskal-Wallis test.
    
    Parameters:
    df (pd.DataFrame): DataFrame containing group labels and AUC values.
    group_column (str): Column name for group labels.
    auc_column (str): Column name for AUC values.
    
    Returns:
    pd.DataFrame: A DataFrame with Dunn's test results.
    """
    # Extract group labels and AUC values
    groups = df[group_column].unique()
    auc_values = [df[df[group_column] == group][auc_column].values for group in groups]
    
    # Perform Kruskal-Wallis test
    _, p_value_kruskal = kruskal(*auc_values)
    #print(f'Kruskal-Wallis test p-value: {p_value_kruskal}')
    
    if p_value_kruskal > 0.05:
        #print("No significant differences between groups (Kruskal-Wallis test).")
        return None
    
    # Perform Dunn's test
    dunns_results = sp.posthoc_dunn(df, val_col=auc_column, group_col=group_column)
    
    return dunns_results


# Get the first column using iloc
first_column = df.iloc[:, 0]

# Count the number of unique groups in the first column
num_groups = first_column.nunique()
group_s = df.shape[0]/num_groups
# Round down to the nearest integer
group_size = math.floor(group_s)


num_iterations = 1000
num_sig = 0
dunns_results_list = []
for i in range(num_iterations):
    # Split the DataFrame into 3 random groups
    groups = split_dataframe_randomly(df, num_groups)
    
    for i, group in enumerate(groups):
        group['AUC'] = calculate_auc_per_row(group)
        
    combined = replace_labels_and_combine(groups, 'Treatment') 
    
    dunns_results = perform_dunns_test(combined, 'Treatment', 'AUC')
    dunns_results_list.append(dunns_results)
    if(isinstance(dunns_results, pd.DataFrame)):
        num_sig = num_sig+1

estimate_fp_rate = num_sig/num_iterations