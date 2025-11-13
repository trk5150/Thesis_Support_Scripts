# Install necessary packages
if (!requireNamespace("rentrez", quietly = TRUE)) install.packages("rentrez")
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")

library(rentrez)
library(Biostrings)

# Function to extract species name from sequence headers
extract_species_from_headers <- function(headers) {
  # Assuming species names are in a consistent format, like "Homo_sapiens" or "Homo sapiens"
  # Adapt this regex as needed for your alignment file
  species_names <- sapply(headers, function(header) {
    match <- regexpr("[A-Z][a-z]+_[a-z]+|[A-Z][a-z]+\\s[a-z]+", header)
    if (match > 0) {
      return(substr(header, match, match + attr(match, "match.length") - 1))
    } else {
      return(NA)
    }
  })
  return(unique(na.omit(species_names))) # Return unique species names
}

# Function to fetch taxonomy using NCBI
fetch_taxonomy <- function(species_name) {
  # Search for the species in the taxonomy database
  search_results <- entrez_search(db = "taxonomy", term = species_name, retmax = 1)
  
  if (length(search_results$ids) > 0) {
    # Fetch taxonomy data using the first ID
    tax_data <- entrez_summary(db = "taxonomy", id = search_results$ids[1])
    taxonomy <- tax_data$lineage
    return(taxonomy)
  } else {
    return(paste("No taxonomy found for", species_name))
  }
}

# Main function to process the FASTA file
process_alignment_file <- function(file_path) {
  # Read the FASTA file
  alignment <- readDNAStringSet(file_path, format = "fasta")
  
  # Extract headers
  headers <- names(alignment)
  
  # Extract species names from headers
  species_names <- extract_species_from_headers(headers)
  
  # Fetch taxonomy for each species
  taxonomy_results <- sapply(species_names, fetch_taxonomy)
  
  # Combine results
  results <- data.frame(Species = species_names, Taxonomy = taxonomy_results)
  
  # Print results
  print(results)
  return(results)
}


# Example usage
file_path <- "C://Users//tik105//OneDrive - Harvard University//Atollin Paper//Conservation//PlotProtein//MUSCLE aligned.fas"  # Path to the FASTA alignment file
taxonomy_results <- process_alignment_file(file_path)
