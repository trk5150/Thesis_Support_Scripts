import pandas as pd
from Bio import pairwise2

# Load the target protein sequence
#target_sequence = "KCPMLRSRLGQESVHCGPMFIQVSRPLPLWRDNRQTPWLLSLRGELVASLEDASLMGLYVDMNATTVTVQSPRQGLLQRWEVLNTSAELLPLWLVSGHHAYSLEAACPPVSFQPESEVLVHIPKQRLGLVKRGSYIEETLSLRFLRVHQSNIFMVTENKDFVVVSIPAAGVLQVQRCQEVGGTPGTQAFYRVDLSLEFAEMAAPVLWTVESFFQC"
target_sequence = "KCPMLRSRLGQESVHCGPMFIQVSRPLPLWRDNRQTPWLLSL"
# Load the Excel file
file_path = "C://Users//tik105//Downloads//neuropeptide_excel_NeuroPepV2_neuropeptide_all.xlsx"
sheet_name = 0  # Adjust if needed

# Read the Excel file
df = pd.read_excel(file_path, sheet_name=sheet_name)

# Define column names
sequence_column = "Sequence"
name_column = "Name"
species_column = "Organism"
size_column = "Length"

# Validate columns
missing_columns = [col for col in [sequence_column, name_column, species_column, size_column] if col not in df.columns]
if missing_columns:
    raise ValueError(f"Missing required columns in Excel file: {', '.join(missing_columns)}")

# Store alignment results
alignment_results = []

# Iterate over sequences and perform alignment
for _, row in df.iterrows():
    protein_seq = str(row[sequence_column]).strip()
    name = str(row[name_column]).strip()
    species = str(row[species_column]).strip()
    size = row[size_column]  # Assuming Length is numeric

    # Skip empty sequences
    if not protein_seq:
        continue

    # Perform pairwise alignment (score-only)
    actual_score = pairwise2.align.globalxx(protein_seq, target_sequence, score_only=True)

    # Compute max possible score by aligning the sequence with itself
    max_possible_score = pairwise2.align.globalxx(protein_seq, protein_seq, score_only=True)

    # Avoid division by zero
    normalized_score = actual_score / max_possible_score if max_possible_score > 0 else 0

    # Store results
    alignment_results.append({
        "Name": name,
        "Organism": species,
        "Length": size,
        "Score": actual_score,
        "Max Score": max_possible_score,
        "Normalized Score": round(normalized_score, 4)  # Rounds to 4 decimal places
    })

# Convert to DataFrame
results_df = pd.DataFrame(alignment_results)

excel_output_file = "C://Users//tik105//Downloads//neuropeptide_alignment_45aa-frag_results.xlsx"
with pd.ExcelWriter(excel_output_file, engine='openpyxl') as writer:
    results_df.to_excel(writer, sheet_name="All Alignments", index = False)
    

