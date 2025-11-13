import pandas as pd
import plotly.express as px
import plotly.io as pio

# Load the Excel file
file_path = "C://Users//tik105//OneDrive - Harvard University//Atollin Paper//Conservation//BlastPResults_Taxonomy - no scores.xlsx"  # Replace with the path to your Excel file
df = pd.read_excel(file_path)

# Define the species to include
defined_species = ["Danio rerio", "Xenopus laevis", "Homo sapiens", "Pan troglodytes", "Macaca mulatta", "Rattus rattus", "Mus musculus"]

# Create a new column for selective species inclusion
df["filtered_species"] = df["Species"].where(df["Species"].isin(defined_species), "")

# Create the sunburst plot
fig = px.sunburst(
    df,
    path=["phylum", "class", "order", "family", "filtered_species"],  # Hierarchical structure
    title="Sunburst Plot with Selective Species Inclusion",
    color="family",  # Color based on phylum
    color_discrete_sequence=px.colors.qualitative.Set3  # Adjust color scheme as needed
)

# Update marker colors to make empty wedges transparent
fig.data[0].marker.colors = [
    "rgba(0,0,0,0)" if species == "" else color
    for species, color in zip(df["filtered_species"], fig.data[0].marker.colors)
]

# Save the sunburst plot as an HTML file
output_path = 'C://Users//tik105//Desktop//Sun_plot.html'  # Specify output path
pio.write_html(fig, output_path)

# Print a confirmation message
print(f"Sunburst plot saved successfully to: {output_path}")
