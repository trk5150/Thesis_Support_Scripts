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
adata = sc.read_h5ad("C://Users//trk51//Downloads//TabSap//Pancreas_TSP1_30_version2d_10X_smartseq_scvi_Nov122024.h5ad")

# View metadata (obs) and features (var)
print(adata.obs.head())  # Cell metadata (samples, cell types, conditions, etc.)
print(adata.var.head())  # Gene metadata

meta2 = adata.obs
meta = adata.var

# Define which column in `obs` contains sample and cell type
sample_col = "ensembl_id"   # Adjust this to your actual column name
celltype_col = "cell_ontology_class"  # Adjust this to your actual column name

# Define the gene of interest
gene_of_interest = "C1orf127"  # Replace with the actual gene name

# Check if the gene exists in the dataset
if gene_of_interest not in adata.var["gene_symbol"].values:
    raise ValueError(f"Gene '{gene_of_interest}' not found in the dataset.")

gene_id = adata.var.index[adata.var["gene_symbol"] == gene_of_interest][0]  # Gets the corresponding gene ID

# Get the index of the gene in the matrix
gene_idx = list(adata.var.index).index(gene_id)  # Converts the index to a list and finds position

# Use sparse slicing to keep memory usage low
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

# Count the number of cells per cell type
cell_counts = gene_df["cell_type"].value_counts().reset_index()
cell_counts.columns = ["cell_type", "cell_count"]

# Merge with pseudobulk gene counts
pseudobulk_gene = pseudobulk_gene.merge(cell_counts, on="cell_type")

# Normalize by the number of cells in each cell type
pseudobulk_gene["normalized_gene_count"] = pseudobulk_gene["gene_count"] / pseudobulk_gene["cell_count"]

# Save to CSV (optional)
pseudobulk_gene.to_csv("C://Users//trk51//Downloads//panc_pseudobulk_normalized_counts.csv", index=False)

