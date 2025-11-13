# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 13:50:26 2025

@author: tik105
"""

import pandas as pd
from Bio import Entrez

# Set your email for Entrez (required by NCBI)
Entrez.email = "tkunz@fas.harvard.edu"

def get_taxonomic_ranks(genus_species, counter):
    print(f"Querying NCBI for: {genus_species}")
    try:
        # Search for taxonomy ID
        handle = Entrez.esearch(db="taxonomy", term=genus_species)
        record = Entrez.read(handle)
        handle.close()
        
        if record["IdList"]:
            tax_id = record["IdList"][0]
            # Fetch taxonomy details
            handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
            record = Entrez.read(handle)
            #print(record)
            handle.close()
            
            # Extract specific ranks
            lineage_ex = record[0].get("LineageEx", [])
            taxonomy = {
                "kingdom": None,
                "phylum": None,
                "class": None,
                "order": None,
                "family": None,
                "genus": None,
                "Species": record[0].get("ScientificName", None),
            }
            for rank in lineage_ex:
                if rank["Rank"] in taxonomy:
                    taxonomy[rank["Rank"]] = rank["ScientificName"]
        else:
            taxonomy = {
                "kingdom": "Not found",
                "phylum": "Not found",
                "class": "Not found",
                "order": "Not found",
                "family": "Not found",
                "genus": "Not found",
                "species": "Not found",
            }
        
        # Print progress for every 25th species
        if counter % 25 == 0:
            print(f"Processed {counter} species so far...")
        
        return taxonomy
    except Exception as e:
        return {
            "Kingdom": f"Error: {e}",
            "Phylum": f"Error: {e}",
            "Class": f"Error: {e}",
            "Order": f"Error: {e}",
            "Family": f"Error: {e}",
            "Genus": f"Error: {e}",
            "Species": f"Error: {e}",
        }

# Load the BLASTP results from Excel
file_path = "C://Users//tik105//OneDrive - Harvard University//Atollin Paper//Conservation//BlastPResults_no_duplicates.xlsx"

df = pd.read_excel(file_path)

# Initialize a counter
counter = 0

# Add new columns for taxonomic ranks
taxonomic_data = []
for index, species in enumerate(df["Scientific Name"], start=1):
    counter += 1
    taxonomy = get_taxonomic_ranks(species, counter)
    taxonomic_data.append(taxonomy)

# Convert the list of taxonomic data into a DataFrame
taxonomy_df = pd.DataFrame(taxonomic_data)

# Merge the new taxonomic columns with the original DataFrame
df = pd.concat([df, taxonomy_df], axis=1)

# Save the updated DataFrame back to Excel
output_file = "C://Users//tik105//OneDrive - Harvard University//Atollin Paper//Conservation//BlastPResults_Taxonomy.xlsx"
df.to_excel(output_file, index=False)

print(f"Updated results saved to {output_file}")
