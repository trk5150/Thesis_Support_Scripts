import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import os

# Define the directory where the .h5ad files are located
data_directory = "C://Users//trk51//Downloads//TabSap//"

# List all .h5ad files in the directory
h5ad_files = [f for f in os.listdir(data_directory) if f.endswith('.h5ad')]

# Initialize an empty list to store the pseudobulk results
pseudobulk_results = []

# Define the column names for the metadata and gene of interest
sample_col = "ensembl_id"   # Adjust this to your actual column name
celltype_col = "tissue"  # Adjust this to your actual column name
#celltype_col = "notes"
gene_of_interest = "C1orf127"  # Change this to the gene you're interested in

count = 0

# Loop over each .h5ad file in the directory
for h5ad_file in h5ad_files:
    # Load the .h5ad file
    adata = sc.read_h5ad(os.path.join(data_directory, h5ad_file))

    # Check if the gene exists in the dataset
    if gene_of_interest not in adata.var["gene_symbol"].values:
        print(f"Gene '{gene_of_interest}' not found in {h5ad_file}. Skipping this file.")
        continue
    
    # Get the gene ID for the gene of interest
    gene_id = adata.var.index[adata.var["gene_symbol"] == gene_of_interest][0]

    # Get the index of the gene in the matrix
    gene_idx = list(adata.var.index).index(gene_id)  

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

    # Count the number of cells per cell type
    cell_counts = gene_df["cell_type"].value_counts().reset_index()
    cell_counts.columns = ["cell_type", "num_cells"]

    # Aggregate total gene counts per cell type
    pseudobulk_gene = gene_df.groupby("cell_type", observed=False)["gene_count"].sum().reset_index()

    # Merge with the cell count data
    pseudobulk_gene = pseudobulk_gene.merge(cell_counts, on="cell_type", how="left")

    # Compute normalized gene expression (average expression per cell)
    pseudobulk_gene["normalized_expression"] = pseudobulk_gene["gene_count"] / pseudobulk_gene["num_cells"]

    # Add the sample file name to track the source
    pseudobulk_gene['source_file'] = h5ad_file

    # Append the results to the list
    pseudobulk_results.append(pseudobulk_gene)
    
    count += 1
    print(f"Processed {count} files.")

# Combine all pseudobulk results into a single DataFrame
combined_pseudobulk = pd.concat(pseudobulk_results, ignore_index=True)

# Save to CSV
combined_pseudobulk.to_csv("C://Users//trk51//Downloads//combined_pseudobulk_counts.csv", index=False)

print("Pseudobulk analysis with normalization completed and saved.")
