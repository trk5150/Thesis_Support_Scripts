import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio

# Step 1: Define the matrix
data_matrix = pd.read_csv("C://Users//tik105//Dropbox (Harvard University)//ER-seq htv//Every GTT collection//Summary analyses//Summary.csv", index_col=0)

# Debug: Check the original DataFrame
print("Original Data Matrix:")
print(data_matrix)

# Step 2: Convert strings to tuples
def str_to_tuple(s):
    return tuple(map(int, s.split(', ')))  # Split by ', ' and convert to integers

# Apply the conversion to the entire DataFrame
data_matrix = data_matrix.applymap(str_to_tuple)

# Step 2: Extract data for Sankey diagram
long_data = []

# Define a color map for each unique label
color_map = {
    'C1orf127 full': '#9D00EB',
    'N-terminal fragment': '#534CFA',
    'C1orf127 fragment other': '#69A2EB',
    'Other constructs': '#808080',
    'HTVI': '#8AE6E1',
    'AAV': '#8ACDE6',
    'Peptide injection': '#446470',
    'Genetic OE': '#8AB0E6',
    'Significant': '#9D00EB',
    'Non-significant': '#808080'
}

# Iterate through the matrix to fill the long_data list
for source in data_matrix.index:
    for target in data_matrix.columns:
        value = data_matrix.at[source, target]  # Access the cell
        if isinstance(value, tuple):  # Ensure the cell contains a tuple
            n, m = value  # Unpack the tuple
            if n > 0 or m > 0:  # Only include non-zero flows
                # Append data for flow from Source to Target (first)
                long_data.append({'source': source, 'target': target, 'value': n + m, 'color': color_map[source]})

                # Append sub-target flows (second)
                long_data.append({'source': target, 'target': "Significant", 'value': n, 'color': color_map[source]})
                long_data.append({'source': target, 'target': "Non-significant", 'value': m, 'color': color_map[source]})

# Convert to DataFrame
long_data_df = pd.DataFrame(long_data)

# Step 3: Create the Sankey diagram
# Get unique labels for nodes
labels = list(pd.concat([long_data_df['source'], long_data_df['target']]).unique())

# Map source and target labels to indices in the labels list
long_data_df['source_index'] = long_data_df['source'].apply(lambda x: labels.index(x))
long_data_df['target_index'] = long_data_df['target'].apply(lambda x: labels.index(x))

# Assign colors to links based on their source (same color for direct and sub-target flows)
link_colors = long_data_df['color'].tolist()

# Step 4: Create the Sankey diagram
fig = go.Figure(go.Sankey(
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color="black", width=0.5),
        label=labels,
        color=[color_map.get(label, '#7f7f7f') for label in labels]  # Use the color map
    ),
    link=dict(
        source=long_data_df['source_index'],
        target=long_data_df['target_index'],
        value=long_data_df['value'],
        color=link_colors  # Assigning colors to links
    )
))

# Step 5: Display the diagram
fig.update_layout(title_text="Sankey Diagram with Custom Color Mapping", font_size=12)
pio.write_html(fig, 'C://Users//tik105//Desktop//GTT false positive estimation//your_plot.html')
