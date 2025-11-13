import pandas as pd

# Load the Excel file
input_file = "C://Users//tik105//Desktop//Diff 9 Drug TC and results//FACs+ELISA combined results.xlsx"
df = pd.read_excel(input_file)

# Clean column names
df.columns = df.columns.str.strip()

# Define group keys and relevant numeric data columns
group_cols = ['Sample or drug', 'Dose']
data_cols = [
    'Cell count', 'Beta count', '% beta',
    'Low glucose', 'High glucose', 'KCl',
    'LG per 1000 cell', 'HG per 1000 cell', 'KCl per cell',
    'LG per 1000 beta', 'HG per 1000 beta', 'KCl per 1000 beta',
    'Stim index'
]

# Group and aggregate
grouped = df.groupby(group_cols)
aggregated = grouped[data_cols].agg(['mean', 'std', 'count'])

# Flatten multi-level columns
aggregated.columns = ['_'.join(col).strip() for col in aggregated.columns.values]
aggregated.reset_index(inplace=True)

# --- NEW: Convert dose to numeric value for sorting ---
def parse_dose(dose):
    if isinstance(dose, (int, float)):
        return float(dose)
    dose = str(dose).lower().strip()
    if dose == '-1':
        return -1.0
    elif dose.endswith('um'):
        return float(dose.replace('um', ''))
    elif dose.endswith('nm'):
        return float(dose.replace('nm', '')) / 1000  # convert nM to uM
    else:
        return float('inf')  # unknown values sorted last

aggregated['Dose_numeric'] = aggregated['Dose'].apply(parse_dose)

# Combine Sample and Dose into one column
aggregated['Sample and Dose'] = aggregated['Sample or drug'].astype(str) + ' ' + aggregated['Dose'].astype(str)

# Write one sheet per data column
output_file = "C://Users//tik105//Desktop//Diff 9 Drug TC and results//aggregated_by_column_combined.xlsx"

with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
    for base_col in data_cols:
        matching_cols = [f"{base_col}_mean", f"{base_col}_std", f"{base_col}_count"]
        selected_cols = ['Sample or drug', 'Dose', 'Dose_numeric', 'Sample and Dose'] + matching_cols
        sheet_data = aggregated[selected_cols].copy()

        # Rename columns for melting
        sheet_data.rename(columns={
            f"{base_col}_mean": "Mean",
            f"{base_col}_std": "Std",
            f"{base_col}_count": "N"
        }, inplace=True)

        # Melt into long format
        long_df = pd.melt(
            sheet_data,
            id_vars=['Sample and Dose', 'Sample or drug', 'Dose', 'Dose_numeric'],
            value_vars=['Mean', 'Std', 'N'],
            var_name='Parameter',
            value_name='Value'
        )

        # Set desired parameter order
        param_order = ['Mean', 'Std', 'N']
        long_df['Parameter'] = pd.Categorical(long_df['Parameter'], categories=param_order, ordered=True)
        
        # Sort using the categorical parameter order
        long_df.sort_values(by=['Sample or drug', 'Dose_numeric', 'Parameter'], inplace=True)


        # Drop sorting helper columns
        long_df.drop(columns=['Sample or drug', 'Dose', 'Dose_numeric'], inplace=True)

        # Write to a safe sheet name
        safe_sheet_name = base_col.replace("/", " per").replace(":", "").strip()[:31]
        long_df.to_excel(writer, sheet_name=safe_sheet_name, index=False)
