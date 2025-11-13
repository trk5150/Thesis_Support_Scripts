import os
import math
import numpy as np
import pandas as pd
from scipy.stats import kruskal
import scikit_posthocs as sp

def split_dataframe_randomly(df, n_groups):
    """Splits a DataFrame into a set number of random groups."""
    n_rows = len(df)
    shuffled_indices = np.random.permutation(n_rows)
    group_sizes = np.array_split(shuffled_indices, n_groups)
    groups = [df.iloc[indices].reset_index(drop=True) for indices in group_sizes]
    return groups

def calculate_auc_per_row(df):
    """Calculates the area under the curve (AUC) for each row."""
    time_points = df.columns[1:].astype(float).values
    time_course_data = df.iloc[:, 1:].values
    auc_values = [np.trapz(y=row, x=time_points) for row in time_course_data]
    return pd.Series(auc_values, index=df.index)

def replace_labels_and_combine(dfs, base_label='Treatment'):
    """Replaces labels and combines DataFrames."""
    combined_df = pd.DataFrame()
    for i, df in enumerate(dfs):
        unique_identifier = f"{base_label}{i+1}"
        df[df.columns[0]] = unique_identifier
        combined_df = pd.concat([combined_df, df], ignore_index=True)
    return combined_df

def perform_dunns_test(df, group_column, auc_column):
    """Performs Dunn's test on AUC values."""
    groups = df[group_column].unique()
    auc_values = [df[df[group_column] == group][auc_column].values for group in groups]
    _, p_value_kruskal = kruskal(*auc_values)
    
    if p_value_kruskal > 0.05:
        return None
    
    dunns_results = sp.posthoc_dunn(df, val_col=auc_column, group_col=group_column)
    return dunns_results

def process_excel_file(file_path):
    """Processes each Excel file and calculates Dunn's test results."""
    df = pd.read_excel(file_path)
    
    first_column = df.iloc[:, 0]
    num_groups = first_column.nunique()
    group_size = math.floor(df.shape[0] / num_groups)
    
    num_iterations = 1000
    num_sig = 0
    dunns_results_list = []
    
    for _ in range(num_iterations):
        groups = split_dataframe_randomly(df, num_groups)
        
        for group in groups:
            group['AUC'] = calculate_auc_per_row(group)
        
        combined = replace_labels_and_combine(groups, 'Treatment') 
        
        dunns_results = perform_dunns_test(combined, 'Treatment', 'AUC')
        dunns_results_list.append(dunns_results)
        
        if isinstance(dunns_results, pd.DataFrame):
            num_sig += 1
    
    estimate_fp_rate = num_sig / num_iterations * 100
    return estimate_fp_rate

def process_all_excel_files(folder_path):
    """Processes all Excel files in the folder."""
    all_results = {}
    
    for file_name in os.listdir(folder_path):
        if file_name.endswith('.xlsx'):
            file_path = os.path.join(folder_path, file_name)
            print(f"Processing file: {file_name}")
            result = process_excel_file(file_path)
            all_results[file_name] = result
    
    return all_results

# Example usage
#folder_path = 'C://Users//tik105//Desktop//GTT false positive estimation/'
folder_path = 'C://Users//tik105//Dropbox//ER-seq htv//Every GTT collection'
results = process_all_excel_files(folder_path)

# Print results
for file_name, estimate_fp_rate in results.items():
    print(f"{estimate_fp_rate:.4f}")
