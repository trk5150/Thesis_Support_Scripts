# Install and load necessary libraries
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install(c("Biostrings", "msa"))
#install.packages(c("ape", "phangorn", "ggtree"))

library(Biostrings)
library(msa)
library(ape)
library(phangorn)
library(ggtree)
library(ggplot2)

# Step 1: Load and Parse Alignment File
alignment_file <- "C://Users//tik105//OneDrive - Harvard University//Atollin Paper//Conservation//Aligned_Model-Orgs_sequence.fas"  # Path to the FASTA alignment file
alignment <- readAAStringSet(alignment_file)  # Read protein alignment as an AAStringSet

# Parse names to extract the last two elements (e.g., "Sapajus apella")
parsed_names <- sapply(names(alignment), function(name) {
  parts <- strsplit(name, " ")[[1]]  # Split name by spaces
  paste(tail(parts, 2), collapse = " ")  # Combine the last two elements
})

# Update the alignment object with parsed names
names(alignment) <- parsed_names

# Step 2: Convert Alignment to AAbin for Distance Calculation
alignment_bin <- as.AAbin(alignment)

# Step 3: Calculate Pairwise Distances
distances <- dist.ml(alignment_bin, model = "JTT")  # Use JTT model for proteins

# Step 4: Build the Initial Phylogenetic Tree
tree <- nj(distances)  # Neighbor-Joining tree based on full alignment

min_positive_length <- min(tree$edge.length[tree$edge.length > 0])
tree$edge.length[tree$edge.length < 0] <- min_positive_length



# Step 6: Convert Tree to Ultrametric Using `chronos`
ultrametric_tree <- chronos(tree, model = "relaxed", lambda = 1)  # Make the tree ultrametric

# Calculate maximum branch length for scaling
max_branch_length <- max(ape::branching.times(ultrametric_tree))

# Plot the tree without the x-axis or the black line
ggtree(ultrametric_tree, layout = "rectangular") +
  geom_tiplab(align = TRUE, linesize = 1, size = 3) +  # Align tip labels and adjust text size
  theme_tree2() +
  xlim(0, max_branch_length + 0.4 * max_branch_length) +  # Add extra space for labels
  theme(
    axis.text.x = element_blank(),        # Remove x-axis text
    axis.ticks.x = element_blank(),       # Remove x-axis ticks
    axis.title.x = element_blank(),       # Remove x-axis title
    panel.grid = element_blank(),         # Remove grid lines
    panel.border = element_blank(),       # Remove the border (black line)
    plot.margin = margin(10, 10, 10, 10) # Add margin for tip labels
  )