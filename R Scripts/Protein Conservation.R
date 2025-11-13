# Install required packages if not already installed
install.packages("BiocManager")
install.packages("plotprotein")


library(BiocManager)

BiocManager::install(update = TRUE, ask = FALSE)

BiocManager::install(c("msa", "seqinr", "Biostrings"))


# Load libraries
library(msa)          # For multiple sequence alignment
library(seqinr)       # For sequence handling
library(Biostrings)   # For working with biological sequences
library(plotprotein)  # For visualization

# Load the FASTA file
sequences <- Biostrings::readAAStringSet("C://Users//tik105//OneDrive - Harvard University//Atollin Paper//Conservation//Non-repeat sequences.fasta")

# Perform MUSCLE alignment
alignment <- msa(sequences, method = "MUSCLE")

# Convert alignment to a matrix for further analysis
alignment_matrix <- msaConvert(alignment, type = "seqinr::alignment")

# Convert alignment to a matrix
aligned_sequences <- as.matrix(alignment_matrix$seq)

# Convert sequences from single strings to individual residues
aligned_sequences <- do.call(rbind, strsplit(as.character(aligned_sequences[, 1]), split = ""))

# Check the dimensions
dim(aligned_sequences)  # Should be rows = sequences, columns = alignment positions

# Calculate entropy-based conservation scores
conservation_scores <- apply(aligned_sequences, 2, function(column) {
  freq <- table(column) / length(column)  # Normalize frequency by column length
  -sum(freq * log2(freq + 1e-9))          # Add a small value to avoid log(0)
})

# Normalize scores (1 = most conserved)
normalized_scores <- 1 - (conservation_scores / max(conservation_scores))

# Prepare conservation data for plotting
conservation_data <- data.frame(
  Position = 1:length(normalized_scores),
  Conservation = normalized_scores
)

ggplot(conservation_data, aes(x = Position, y = Conservation)) +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Conservation Trace",
    x = "Residue Position",
    y = "Conservation Score"
  )
# Example: Define domains
domains <- data.frame(
  start = c(10, 50, 100),
  end = c(30, 70, 120),
  label = c("Domain1", "Domain2", "Domain3")
)


## based on similarity scores

# Extract the reference sequence
reference_sequence <- aligned_sequences[1, ]

# Calculate similarity scores
similarity_scores <- apply(aligned_sequences, 2, function(column) {
  sum(column == reference_sequence) / length(column)  # Fraction of sequences matching the reference
})
normalized_similarity <- similarity_scores / max(similarity_scores)

conservation_data <- data.frame(
  Position = 1:length(normalized_similarity),
  Similarity = normalized_similarity
)

library(ggplot2)

ggplot(conservation_data, aes(x = Position, y = Similarity)) +
  geom_line(color = "blue") +
  theme_minimal() +
  labs(
    title = "Conservation Relative to Reference Sequence",
    x = "Residue Position",
    y = "Similarity Score"
  )



