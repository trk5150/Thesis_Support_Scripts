library(readxl)
library(DESeq2)
library(ggplot2) 
library(ggrepel)
library(pheatmap)

# Load your data
data <- read_excel("C://Users//tik105//Desktop//Neonate phenotyping 5-1-2024//Bulk RNA-seq//Null Livers//Results//gene_readcounts_of_all_samples.xlsx")

# Assuming "Gene_Name" column contains gene identifiers
count_data <- as.matrix(data[, c("7-1", "7-2", "7-5", "7-6", "13-1", "13-2", "13-5", "13-7")])
rownames(count_data) <- data$Gene_Name

# Create coldata with batch information
coldata <- data.frame(
  condition = c("Het", "Het", "Null", "Null", "Het", "Null", "Het", "Null"),
  batch = c("7", "7", "7", "7", "13", "13", "13", "13"),  # Batch variable
  row.names = colnames(count_data)
)
coldata$condition <- factor(coldata$condition)  # Ensure condition is a factor
coldata$batch <- factor(coldata$batch)          # Ensure batch is a factor

# Create DESeqDataSet and include batch in the design
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = coldata, design = ~ batch + condition)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results for the condition effect (adjusted for batch)
res <- results(dds, contrast = c("condition", "Null", "Het"))

# Sort results by adjusted p-value
res <- res[order(res$padj), ]

# View the top results
head(res)

# Save results to a CSV file
write.csv(as.data.frame(res), "C://Users//tik105//Desktop//Neonate phenotyping 5-1-2024//Bulk RNA-seq//Null Livers//Results//DiffEx_with_batch_correction.csv")

# MA Plot
plotMA(res, main = "DESeq2 (Batch Corrected)", ylim = c(-2, 2))

# Volcano Plot
volcano_data <- as.data.frame(res)
volcano_data$significance <- ifelse(volcano_data$padj < 0.05 & volcano_data$log2FoldChange > 0, "Upregulated",
                                    ifelse(volcano_data$padj < 0.05 & volcano_data$log2FoldChange < 0, "Downregulated", "Non-significant"))

# Specify your genes of interest
genes_of_interest <- c("G6pc")  # Replace with your genes of interest
# Specify genes of interest: top 20 genes for both positive and negative fold changes
top_positive <- rownames(volcano_data[volcano_data$log2FoldChange > 0, ][order(volcano_data$padj[volcano_data$log2FoldChange > 0]), ])[1:10]
top_negative <- rownames(volcano_data[volcano_data$log2FoldChange < 0, ][order(volcano_data$padj[volcano_data$log2FoldChange < 0]), ])[1:10]
# Combine top positive and negative genes along with any genes of interest
all_genes_to_label <- unique(c(top_positive, top_negative, genes_of_interest))
# Add labels for the selected genes
volcano_data$label <- ifelse(rownames(volcano_data) %in% all_genes_to_label, rownames(volcano_data), NA)

# Create the volcano plot
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), size = 1.5) +
  scale_color_manual(values = c("Non-significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) +
  xlim(-5, 5) +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value", title = "Volcano Plot of Differential Expression (Batch Corrected)") +
  theme_minimal() +
  theme(legend.position = "top") +
  geom_text_repel(aes(label = label), size = 3, color = "black", max.overlaps = 50, na.rm = TRUE)

# Heatmap of top 100 significant genes
top_genes <- rownames(res[order(res$padj), ])[1:100]
normalized_counts <- counts(dds, normalized = TRUE)
top_counts <- normalized_counts[top_genes, ]

pheatmap(
  top_counts,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  annotation_col = data.frame(condition = coldata$condition),
  scale = "row",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Heatmap of Top 100 Genes (Batch Corrected)",
  treeheight_row = 10,
  treeheight_col = 10,
  border_color = NA,
  fontsize = 8,
  fontsize_row = 6,
  fontsize_col = 8
)
