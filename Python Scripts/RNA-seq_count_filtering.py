import pandas as pd

# Load the Excel file into a DataFrame
def read_excel_file(excel_path):
    df = pd.read_excel(excel_path)
    return df

# Read the gene names from the line-split list file
def read_gene_names(file_path):
    with open(file_path, 'r') as file:
        gene_names = file.read().splitlines()  # Read and split by line
    return gene_names

# Filter the DataFrame to include only rows with matching Gene_Name in the gene list
def filter_dataframe_by_gene_list(df, gene_names):
    # Convert the 'Gene_Name' column and the gene_names list to lowercase
    filtered_df = df[df['Gene_Name'].str.lower().isin([name.lower() for name in gene_names])]
    return filtered_df

# Main function to handle the process
def process_files(excel_file_path, gene_list_file_path, output_file_path):
    df = read_excel_file(excel_file_path)
    gene_names = read_gene_names(gene_list_file_path)
    filtered_df = filter_dataframe_by_gene_list(df, gene_names)

    # Save the filtered DataFrame to a new Excel file
    filtered_df.to_excel(output_file_path, index=False)

# Example usage:
excel_file = "D://Bulk RNA-seq results//24231-02-Analysis-09262024_142004 (neonate livers)//Base counts.xlsx"           # Path to your Excel file
gene_list_file =   "D://Bulk RNA-seq results//24231-02-Analysis-09262024_142004 (neonate livers)//HP_CIRRHOSIS.v2024.1.Hs.grp"    # Path to your gene list file
output_file = "D://Bulk RNA-seq results//24231-02-Analysis-09262024_142004 (neonate livers)//Cirrhosis.xlsx" # Output file path for the filtered data

process_files(excel_file, gene_list_file, output_file)
