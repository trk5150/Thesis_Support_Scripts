# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 17:11:03 2024

@author: tik105
"""


import os
from collections import Counter
import csv

# Directory where your files are stored
directory = "C://Users//tik105//Desktop//LINCs sets//GMT parsed only small molecules"

# Initialize a Counter to store the occurrences of each string
string_counter = Counter()

# Loop through the files in the directory
for filename in os.listdir(directory):
    if filename.endswith(".txt"):  # Only process .txt files
        # Split the filename by '_'
        parts = filename.split('_')
        if len(parts) > 2:  # Ensure there are enough parts in the filename
            string_between = parts[1]  # Extract the string between first and second '_'
            string_counter[string_between] += 1  # Increment the count for this string

# Define the path for the CSV file
csv_file_path = "C://Users//tik105//Desktop//LINCs sets//string_counts.csv"

# Save the counter to a CSV file
with open(csv_file_path, mode='w', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)
    
    # Write the header row
    csv_writer.writerow(['String Between First and Second "_"', 'Occurrences'])
    
    # Write the counts
    for string, count in string_counter.items():
        csv_writer.writerow([string, count])

print(f"String counts have been saved to {csv_file_path}")
