# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 16:40:28 2025

@author: tik105
"""

import pandas as pd
from graphviz import Digraph

# Load the Excel file
file_path = "C://Users//tik105//OneDrive - Harvard University//Atollin Paper//Conservation//BlastPResults_Taxonomy - no scores.xlsx"  # Replace with your file path
df = pd.read_excel(file_path)

# Initialize the graph
flowchart = Digraph(format='png', strict=True)
flowchart.attr(rankdir='LR', size='8,5')

# Add edges to the graph
for _, row in df.iterrows():
    hierarchy = row[["kingdom", "phylum", "class", "order", "family", "genus", "Species"]]
    for i in range(len(hierarchy) - 1):
        parent, child = hierarchy[i], hierarchy[i + 1]
        if pd.notna(parent) and pd.notna(child):  # Ensure no missing values
            flowchart.edge(parent, child)

# Save and visualize the flowchart
output_file = output_file = "C://Users//tik105//Desktop//taxonomy_flowchart"
flowchart.render(output_file, view=True)
