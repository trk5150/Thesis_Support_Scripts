# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 15:39:00 2024

@author: tik105
"""

import pandas as pd
from scipy.stats import hypergeom

# Load the first file (list of chemical compounds by name)
compounds_df = pd.read_csv("C://Users//tik105//OneDrive - Harvard University//SOMSCAN draft//SOMSCAN figures//Drug type overrepressentation testing//Hits list no duplicates.txt")  # adjust the path and file format if necessary

# Load the second file (list of compounds with genes targeted and other data)
targets_df = pd.read_excel("C://Users//tik105//OneDrive - Harvard University//SOMSCAN draft//SOMSCAN figures//Drug type overrepressentation testing//drug metadata.xlsx")

# Assuming both DataFrames have a common column named 'Compound' for the compound names
# If the columns have different names, adjust accordingly
merged_df = pd.merge(compounds_df, targets_df[['Compound', 'Target']], on='Compound', how='left')

# Split the 'Target' column by commas and expand into separate rows
merged_df['Target'] = merged_df['Target'].str.split(', ')
exploded_df = merged_df.explode('Target')

targets_df['Target'] = targets_df['Target'].str.split(', ')
exploded_targets = targets_df.explode('Target')

# Count unique targets in 'exploded_df' (merged file)
merged_target_counts = exploded_df['Target'].value_counts()

# Count unique targets in 'exploded_targets' (whole targets file)
total_target_counts = exploded_targets['Target'].value_counts()

# Combine the counts into a single DataFrame for comparison
comparison_df = pd.DataFrame({
    'Merged_Target_Counts': merged_target_counts,
    'Total_Target_Counts': total_target_counts
}).fillna(0)  # fill NaN values with 0 where targets are missing in one of the files

# Reset index for easier handling and save to a CSV file if needed
comparison_df = comparison_df.reset_index().rename(columns={'index': 'Target'})

target_sig_counts = exploded_targets.groupby('Target')['sig_count'].sum().reset_index()
target_sig_counts = target_sig_counts.rename(columns={'sig_count': 'Total_Sig_Count'})
comparison_df = comparison_df.merge(target_sig_counts, on='Target', how='left')

##Hypergeometric testing:
# Define total population and sample sizes
N = exploded_targets['sig_count'].sum()  # Total number of tests (sum of sig_count)
n = exploded_df['Target'].count()        # Total targets in the merged list (sample size)

# Perform the hypergeometric test for each target
p_values = []
for _, row in comparison_df.iterrows():
    target = row['Target']
    M = row['Total_Sig_Count']           # Total occurrences of the target weighted by sig_count
    k = row['Merged_Target_Counts']      # Observed occurrences of the target in the merged list
    
    # Calculate p-value using the hypergeometric test
    p_value = hypergeom.sf(k - 1, N, M, n)  # sf (survival function) gives P(X >= k)
    p_values.append(p_value)

# Add the p-values to the DataFrame
comparison_df['Hypergeometric_p_value_sigs'] = p_values

# Perform the hypergeometric test for each target
p_values = []
for _, row in comparison_df.iterrows():
    target = row['Target']
    M = row['Total_Target_Counts']           # Total occurrences of the target weighted by sig_count
    k = row['Merged_Target_Counts']      # Observed occurrences of the target in the merged list
    
    # Calculate p-value using the hypergeometric test
    p_value = hypergeom.sf(k - 1, N, M, n)  # sf (survival function) gives P(X >= k)
    p_values.append(p_value)

# Add the p-values to the DataFrame
comparison_df['Hypergeometric_p_value_targets'] = p_values

# Group exploded_targets by 'Target' to get a list of compounds associated with each target
compound_lists = (
    exploded_df
    .groupby('Target')['Compound']
    .unique()
    .apply(lambda x: ', '.join(map(str, x)))  # Convert each item in the list to a string before joining
    .reset_index()
)

compound_lists = compound_lists.rename(columns={'Compound': 'Associated_Compounds'})

# Merge this list with the comparison DataFrame
comparison_df = comparison_df.merge(compound_lists, on='Target', how='left')




# Save the updated DataFrame with the target information
comparison_df.to_csv('C://Users//tik105//OneDrive - Harvard University//SOMSCAN draft//SOMSCAN figures//Drug type overrepressentation testing//hypergeom results no duplicates.csv', index=False)
