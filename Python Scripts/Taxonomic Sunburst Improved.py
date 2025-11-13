import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
# Define the species to include
included_species = [
    "Danio rerio", "Xenopus laevis", "Homo sapiens",
    "Pan troglodytes", "Macaca mulatta", "Rattus rattus", "Mus musculus"
]

# Load data from Excel
file_path = "C://Users//tik105//OneDrive - Harvard University//Atollin Paper//Conservation//BlastPResults_Taxonomy - no scores.xlsx"  # Replace with your file path
data = pd.read_excel(file_path)

# Debugging: Check columns in the loaded file
print("Columns in the Excel file:", list(data.columns))

# Ensure only the relevant columns are selected
expected_columns = ["phylum", "class", "order", "family", "Species"]
relevant_columns = [col for col in expected_columns if col in data.columns]

# Keep only the relevant columns
data = data[relevant_columns]

# Check if all levels are present
if len(relevant_columns) < len(expected_columns):
    print("Warning: Some taxonomic levels are missing in the dataset!")

# Prepare labels, parents, and values for the sunburst chart
sunburst_data = []

# Loop through each level in the hierarchy
for i, level in enumerate(relevant_columns):
    unique_values = data[level].unique()
    for value in unique_values:
        if pd.notna(value):  # Exclude empty or NaN values
            # Identify the parent level
            parent = (
                data.loc[data[level] == value, relevant_columns[i - 1]].iloc[0]
                if i != 0 else ""
            )
            # Aggregate all species under the current taxonomic level
            total_value = data.loc[data[level] == value].shape[0]

            # Add the data to the sunburst_data list
            sunburst_data.append({"label": value, "parent": parent, "value": total_value})

# Convert to DataFrame
sunburst_df = pd.DataFrame(sunburst_data)

# Adjust species values: Set species values to 0 if not in the included list
sunburst_df.loc[
    (sunburst_df["label"].isin(data["Species"])) & (~sunburst_df["label"].isin(included_species)),
    "value"
] = 0

# Verify the Sunburst DataFrame
print("Sunburst DataFrame:")
print(sunburst_df)

# Create the sunburst chart
fig = go.Figure(
    go.Sunburst(
        labels=sunburst_df["label"],
        parents=sunburst_df["parent"],
        values=sunburst_df["value"],
        branchvalues="total"  # Propagate values across the hierarchy
    )
)

import numpy as np

# Apply a logarithmic transformation to the value column
sunburst_df["log_value"] = np.log1p(sunburst_df["value"])  # log1p handles log(0) gracefully

# Update the sunburst chart with the transformed values for coloring
fig.update_traces(
    marker=dict(
        colors=sunburst_df["log_value"],  # Use the transformed 'log_value' column for coloring
        colorscale="Purples",  # Choose a colorscale
        colorbar=dict(title="Log Value")  # Optional: Add a colorbar with a title
    )
)




# Update layout for better display
fig.update_layout(margin=dict(t=0, l=0, r=0, b=0))

# Save the sunburst plot as an HTML file
output_path = 'C://Users//tik105//Desktop//Sun_plot.html'  # Specify output path
pio.write_html(fig, output_path)