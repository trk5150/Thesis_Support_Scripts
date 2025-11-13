# -*- coding: utf-8 -*-

import math
import numpy as np
import pandas as pd
from scipy.stats import kruskal
from scipy.stats import ttest_ind
import scikit_posthocs as sp

# Replace 'file_path.xlsx' with your actual file path
#file_path = "C://Users//tik105//Dropbox (Harvard University)//ER-seq htv//Every GTT collection//2018-3-23 HTVI TK pos.xlsx"
#file_path = "C://Users//tik105//Dropbox (Harvard University)//ER-seq htv//Every GTT collection//2018-5-24 HTVI TK pos.xlsx"
file_path = "C://Users//tik105//Dropbox (Harvard University)//ER-seq htv//Every GTT collection//2024-9-19, genetic OE MALE, IPGTT, TK.xlsx"
df = pd.read_excel(file_path)

def t0_correct(df):
    
    time_columns = df.columns.difference(['Treatment'])

    # Subtract column '0' from all other time columns for each row
    [time_columns] = df[time_columns].sub(df['0'], axis=0)
    df[time_columns] = df[time_columns].sub(df['0'], axis=0)

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
    #print(time_points)
    
    # Extract the time course data (all except the first column, which is 'Treatment')
    time_course_data = df.iloc[:, 1:].values
    #print(time_course_data)
    
    # Calculate AUC for each row using the trapezoidal rule
    auc_values = [np.trapz(y=row, x=time_points) for row in time_course_data]
    
    # Return the AUC values as a pandas Series
    return pd.Series(auc_values, index=df.index)


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
    print(f'Kruskal-Wallis test p-value: {p_value_kruskal}')
    
    #if p_value_kruskal > 0.05:
     #   print("No significant differences between groups (Kruskal-Wallis test).")
      #  return None
    
    # Perform Dunn's test
    dunns_results = sp.posthoc_dunn(df, val_col=auc_column, group_col=group_column)
    print(dunns_results)
    return dunns_results
#def perform_t_test(df, group_column, AuC_column):
    

# Get the first column using iloc
first_column = df.iloc[:, 0]

# Count the number of unique groups in the first column
num_groups = first_column.nunique()
group_s = df.shape[0]/num_groups
# Round down to the nearest integer
group_size = math.floor(group_s)


dunns_results_list = []
df.columns = df.columns.map(str)

df['AUC'] = calculate_auc_per_row(df)

if num_groups>2:
    #peform Tukey's multiple testing if groups >2
    dunns_results = perform_dunns_test(df, 'Treatment', 'AUC')
    
else:
    #perform unpaired t-test if groups<2
    unique_treatments = df['Treatment'].unique()
    # Ensure there are exactly two groups for comparison
    if len(unique_treatments) != 2:
        raise ValueError("There should be exactly two unique treatment groups for an unpaired t-test.")
    
    # Splitting the DataFrame based on the unique treatment groups
    group_1 = df[df['Treatment'] == unique_treatments[0]]
    group_2 = df[df['Treatment'] == unique_treatments[1]]
    
    # List of columns to perform the t-test on (excluding 'Treatment')
    variables = df.columns.drop('Treatment')
    
    # Perform t-tests for each variable
    for var in variables:
        t_stat, p_value = ttest_ind(group_1[var], group_2[var], equal_var=True)
        print(f"T-test for {var}: t-stat = {t_stat}, p-value = {p_value}")
    