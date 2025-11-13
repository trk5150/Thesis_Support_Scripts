library(readxl)
library(DESeq2)
library(ggplot2) 
library(ggrepel)
library(pheatmap)

# Load your data
#data <- read_excel("D://Bulk RNA-seq results//24231-02-Analysis-09262024_142004 (neonate livers)//Base counts.xlsx")
#data <- read_excel("C://Users//tik105//Desktop//Neonate phenotyping 5-1-2024//Bulk RNA-seq//Null Livers//Results//gene_readcounts_of_all_samples.xlsx")
data <- read_excel("C://Users//tik105//Desktop//Neonate phenotyping 5-1-2024//Bulk RNA-seq//Null pancreas//Results//Counts.xlsx")
# Assuming "Gene_Name" column contains gene identifiers
count_data <- as.matrix(data[, c("7-1", "7-2", "7-5", "13-1", "13-2", "13-5", "13-7")])
rownames(count_data) <- data$Gene_Name

coldata <- data.frame(
  condition = c("Het", "Het", "Null", "Het", "Null", "Het", "Null"),
  row.names = colnames(count_data)
)
coldata$condition <- factor(coldata$condition)  # Ensure it's a factor

dds <- DESeqDataSetFromMatrix(countData = count_data, colData = coldata, design = ~ condition)
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "Null", "Het"))

# Sort results by adjusted p-value
res <- res[order(res$padj), ]

# Filter out rows with NA values in the results
res <- na.omit(res)


# View the top results
head(res)

# Save results to a CSV file
write.csv(as.data.frame(res), "C://Users//tik105//Desktop//Neonate phenotyping 5-1-2024//Bulk RNA-seq//Null pancreas//Results//DiffEx.csv")

# MA Plot
plotMA(res, main="DESeq2", ylim=c(-2,2))


##for volcano plot
# Prepare the data for plotting
volcano_data <- as.data.frame(res)
volcano_data$significance <- ifelse(volcano_data$padj < 0.05 & volcano_data$log2FoldChange > 0, "Upregulated",
                                    ifelse(volcano_data$padj < 0.05 & volcano_data$log2FoldChange < 0, "Downregulated", "Non-significant"))

# Specify your genes of interest
genes_of_interest <- c("GM572")  # Replace with your genes of interest
# Specify genes of interest: top 20 genes for both positive and negative fold changes
top_positive <- rownames(volcano_data[volcano_data$log2FoldChange > 0, ][order(volcano_data$padj[volcano_data$log2FoldChange > 0]), ])[1:25]
top_negative <- rownames(volcano_data[volcano_data$log2FoldChange < 0, ][order(volcano_data$padj[volcano_data$log2FoldChange < 0]), ])[1:25]
# Combine top positive and negative genes along with any genes of interest
all_genes_to_label <- unique(c(top_positive, top_negative, genes_of_interest))
# Add labels for the selected genes
volcano_data$label <- ifelse(rownames(volcano_data) %in% all_genes_to_label, rownames(volcano_data), NA)

# Create the volcano plot with customizations
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), size = 1.5) +
  scale_color_manual(values = c("Non-significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) +
  xlim(-5, 5) +  # Limit x-axis to -5 to 5
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value", title = "Volcano Plot of Differential Expression") +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_text_repel(aes(label = label), size = 3, color = "black", max.overlaps = 50, na.rm = TRUE)  # Adjust max.overlaps as needed

##for heatmap:

# Filter significant genes
# Filter the top 50 upregulated and top 50 downregulated genes by padj
top_genes_up <- rownames(res[res$log2FoldChange > 0, ][order(res$padj[res$log2FoldChange > 0]), ])[1:50]
top_genes_down <- rownames(res[res$log2FoldChange < 0, ][order(res$padj[res$log2FoldChange < 0]), ])[1:50]

# Combine the two lists of genes
top_genes <- unique(c(top_genes_up, top_genes_down))

# Subset the results for these top genes
top_genes_res <- res[top_genes, ]

# Extract normalized counts for the selected genes
normalized_counts <- counts(dds, normalized = TRUE)
top_counts <- normalized_counts[top_genes, ]




# Reorder the columns to group Het and Null samples
column_order <- c("7-1", "7-2", "13-5", "13-1",  # Het columns
                  "7-5", "13-2", "13-7")  # Null columns
top_counts <- top_counts[, column_order]  # Reorder columns in the counts matrix

# Update the column annotation to match the reordered columns
annotation_col <- data.frame(
  condition = colData(dds)[column_order, "condition"],
  row.names = column_order
)

# Generate the heatmap with reordered columns
pheatmap(
  top_counts,
  cluster_rows = TRUE,        # Cluster genes by similarity
  cluster_cols = FALSE,       # Do not cluster samples
  show_rownames = TRUE,       # Display gene names as row labels
  annotation_col = annotation_col,  # Annotate columns with condition
  scale = "row",              # Z-score standardization by row (gene)
  color = colorRampPalette(c("blue", "white", "red"))(50),  # Color scheme
  main = "Heatmap of Top 100 Significantly Differentially Expressed Genes",
  treeheight_row = 8,         # Adjust row dendrogram height
  treeheight_col = 8,         # Adjust column dendrogram width
  border_color = NA,          # Remove cell borders for a cleaner look
  fontsize = 7,               # Set base font size for all elements
  fontsize_row = 4,           # Reduce font size for row labels
  fontsize_col = 8            # Set font size for column labels
)

