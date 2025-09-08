import os
import pandas as pd
import re

def parse_subdirectory_column(df):
    """Parses the 'Subdirectory' column into 'Litter', 'Pup', 'Slide', and 'Image Number'."""
    
    litter_list = []
    pup_list = []
    slide_list = []
    image_number_list = []

    for subdir in df['Subdirectory']:
        # Use regex to extract Litter, Pup, Slide, and Image Number
        match = re.search(r'Litter (\d+) pup (\d+) slide (\d+)-(\d+)', subdir)

        if match:
            litter_list.append(int(match.group(1)))  # Extract Litter number
            pup_list.append(int(match.group(2)))     # Extract Pup number
            slide_list.append(int(match.group(3)))   # Extract Slide number
            image_number_list.append(int(match.group(4)))  # Extract first integer after the first dash
        else:
            # If the pattern doesn't match, fill with NaN
            litter_list.append(None)
            pup_list.append(None)
            slide_list.append(None)
            image_number_list.append(None)

    # Add extracted values as new columns
    df['Litter'] = litter_list
    df['Pup'] = pup_list
    df['Slide'] = slide_list
    df['Image Number'] = image_number_list

    return df

def clean_dataframe(df):
    """Removes rows where 'Filename' contains 'c1-4'."""
    
    # Remove rows where 'Filename' column contains 'c1-4'
    df_cleaned = df[~df['Filename'].str.contains('c1-4', na=False)]

    return df_cleaned

def assign_signal(df):
    """Assigns 'Signal' values based on the 'Channel' column."""
    
    signal_mapping = {1: "insulin", 2: "amylase", 3: "glucagon"}
    df['Signal'] = df['Channel'].map(signal_mapping)

    return df

def process_excel_file(excel_path):
    """Loads the Excel file, processes and cleans the DataFrame, assigns signals, and saves the updated file."""
    
    # Load the Excel file
    df = pd.read_excel(excel_path)

    # Process the 'Subdirectory' column
    df = parse_subdirectory_column(df)

    # Clean the DataFrame by removing unwanted rows
    df = clean_dataframe(df)

    # Assign Signal values based on Channel
    df = assign_signal(df)

    # Save the updated DataFrame back to Excel
    updated_excel_path = excel_path.replace(".xlsx", "_cleaned_signal.xlsx")
    df.to_excel(updated_excel_path, index=False)

    print(f"Updated Excel file saved to: {updated_excel_path}")

    return df

# Example usage
excel_path = "D://amylas, ins and gcg costains/Final_Analysis_Results.xlsx"
df_updated = process_excel_file(excel_path)

