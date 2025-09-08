# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 16:57:26 2025

@author: tik105
"""

import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
# Load your single-cell data
adata = sc.read_h5ad("C://Users//tik105//Downloads//2a37a272-6f79-436b-ae11-c1b0b1f293dd.h5ad")

# View metadata (obs) and features (var)
print(adata.obs.head())  # Cell metadata (samples, cell types, conditions, etc.)
print(adata.var.head())  # Gene metadata

meta2 = adata.obs
meta = adata.var

# Define which column in `obs` contains sample and cell type
sample_col = "observation_joinid"   # Change this to match your metadata
#celltype_col = "cell_type" # Change this to match your metadata
celltype_col = "tissue_FACS_droplet" # Change this to match your metadata

# Define the gene of interest
gene_of_interest = "Gm572"  # Replace with the actual gene name

# Check if the gene exists in the dataset
if gene_of_interest not in adata.var["gene_symbols"].values:
    raise ValueError(f"Gene '{gene_of_interest}' not found in the dataset.")

gene_id = adata.var.index[adata.var["gene_symbols"] == gene_of_interest][0]  # Gets the corresponding gene ID

# Get the index of the gene in the matrix
gene_idx = list(adata.var.index).index(gene_id)  # Converts the index to a list and finds position


# Use sparse slicing to keep memory usage low
gene_counts = adata.X[:, gene_idx]

# Ensure sparse matrix slicing works correctly
if sp.issparse(adata.X):
    gene_counts = adata.X[:, gene_idx].toarray().flatten()
else:
    gene_counts = adata.X[:, gene_idx]

# Create a DataFrame with cell types and their corresponding gene counts
gene_df = pd.DataFrame({
    "cell_type": adata.obs[celltype_col].values,
    "gene_count": gene_counts
})

# Aggregate counts at the cell type level
pseudobulk_gene = gene_df.groupby("cell_type")["gene_count"].sum().reset_index()

# Save to CSV (optional)
pseudobulk_gene.to_csv("C://Users//tik105//Downloads//pseudobulk_counts.csv", index=False)
