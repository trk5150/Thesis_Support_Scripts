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

# Step 4: Build the Phylogenetic Tree
tree <- nj(distances)  # Neighbor-Joining tree based on full alignment

# Step 5: Visualize the Tree
# Update tree tip labels with parsed names
tree$tip.label <- names(alignment)


min_positive_length <- min(tree$edge.length[tree$edge.length > 0])
tree$edge.length[tree$edge.length < 0] <- min_positive_length
tree <- ladderize(tree, right = TRUE)
plot.phylo(tree, type = "phylogram", show.tip.label = TRUE)



# Basic visualization
#plot.phylo(tree, type = "phylogram", show.tip.label = TRUE)

# Advanced visualization with ggtree
ggtree(tree) +
  geom_tiplab() +
  theme_tree2()
