import pandas as pd

# Load the Excel file (Update with the actual filename and path)
file_path = "D://amylas, ins and gcg costains//Final_Analysis_Results_cleaned_signal.xlsx"  # Change this to the correct file path
df = pd.read_excel(file_path)

# Drop unnecessary columns
df = df.drop(columns=["Subdirectory", "Filename", "Channel", "Threshold", "Mask File"])

# Pivot the table to collapse rows by Litter, Pup, Slide, and Image Number
df_pivot = df.pivot_table(
    index=["Litter", "Pup", "Slide", "Image Number"],
    columns="Signal",
    values="Stained Area",
    aggfunc="first"
).reset_index()

# Rename columns for clarity
df_pivot.columns.name = None  # Remove hierarchical column name
df_pivot = df_pivot.rename(columns={
    "Amylase": "Stained_Area_Amylase", 
    "Glucagon": "Stained_Area_Glucagon", 
    "Insulin": "Stained_Area_Insulin"
})

# Calculate ratios
df_pivot["Ratio_Insulin_Amylase"] = df_pivot["insulin"] / df_pivot["amylase"]
df_pivot["Ratio_Glucagon_Amylase"] = df_pivot["glucagon"] / df_pivot["amylase"]

# Handle cases where Amylase value is zero to avoid division errors
df_pivot.replace([float('inf'), -float('inf')], None, inplace=True)

# Save the processed DataFrame to a new Excel file
output_file = "D://amylas, ins and gcg costains//collapsed_data.xlsx"
df_pivot.to_excel(output_file, index=False)

print(f"Processed data saved to {output_file}")
