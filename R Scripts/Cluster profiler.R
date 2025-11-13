install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("DOSE")

# Install BiocManager if it's not installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor packages
#BiocManager::install()

library(org.Mm.eg.db)
library(DOSE)
library(clusterProfiler)
library(readxl)
library(DESeq2)
library(ggplot2) 
library(ggrepel)
library(pheatmap)
library(dplyr)
library(stringr)

# Load your data
#data <- read_excel("D://Bulk RNA-seq results//24231-02-Analysis-09262024_142004 (neonate livers)//Base counts.xlsx")
data <- read_excel("C://Users//tik105//Desktop//Neonate phenotyping 5-1-2024//Bulk RNA-seq//Null Livers//Results//gene_readcounts_of_all_samples.xlsx")


#data <- read_excel("D://Bulk RNA-seq results//127-Cag-rosa x ins-cre RNA-seq results//Fat.xlsx")
# Assuming "NAME" column contains gene identifiers
count_data <- as.matrix(data[, c("7-1", "7-2", "7-5", "7-6", "13-1", "13-2", "13-5", "13-7")])
#rownames(count_data) <- data$NAME
rownames(count_data) <- data$Gene_Name

coldata <- data.frame(
  condition = c("Het", "Het", "Null", "Null","Het", "Null", "Het", "Null"),
  row.names = colnames(count_data)
)
coldata$condition <- factor(coldata$condition)  # Ensure it's a factor

dds <- DESeqDataSetFromMatrix(countData = count_data, colData = coldata, design = ~ condition)
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "Het", "Null"))
# Filter out rows with NA values in the results
res <- na.omit(res)

# Sort results by adjusted p-value
res <- res[order(res$padj), ]

# Extract fold changes and store them in a named vector
fold_changes <- res$log2FoldChange
names(fold_changes) <- rownames(res)

# Define thresholds for significance
alpha <- 0.001  # Adjusted p-value threshold (0.05)
lfc_threshold <- 0.5  # Log2 fold change threshold (0.5 = 2x expression difference)
# Remove rows with NA in padj column
res_filtered <- res[!is.na(res$padj), ]

# Extract upregulated genes (Fold change > 1 and padj < alpha)
upregulated_genes <- rownames(res_filtered[res_filtered$padj < alpha & res_filtered$log2FoldChange > lfc_threshold, ])

# Extract downregulated genes (Fold change < -1 and padj < alpha)
downregulated_genes <- rownames(res_filtered[res_filtered$padj < alpha & res_filtered$log2FoldChange < -lfc_threshold, ])

# Perform GO enrichment analysis
go_results <- enrichGO(gene = upregulated_genes,
                       OrgDb = org.Mm.eg.db,  # Change to appropriate organism database
                       keyType = "SYMBOL",   # Change if using different ID type
                       ont = "BP",             # Can be "BP", "CC", or "MF"
                       pAdjustMethod = "BH",   # p-value adjustment method
                       qvalueCutoff = 0.05,    # q-value cutoff
                       readable = TRUE)         # Convert to gene symbols



#head(go_results)
#cnetplot(go_results, foldChange = fold_changes[upregulated_genes])
#barplot(go_results, showCategory = 10)


go_results_df <- as.data.frame(go_results)

# Limit to top 10 results based on Count
top_go_results <- go_results_df %>%
  slice_max(order_by = Count, n = 10) %>%
  mutate(
    Description = ifelse(
      nchar(Description) > 70, 
      str_sub(Description, 1, str_locate(Description, ",")[,1] - 1), 
      Description
    )
  )

# Create new columns for Numerator, Denominator, and GeneRatio Result
top_go_results <- top_go_results %>%
  mutate(
    Numerator = as.numeric(str_extract(GeneRatio, "^[^/]+")),  # Extract numerator
    Denominator = as.numeric(str_extract(GeneRatio, "(?<=/)[0-9]+")),  # Extract denominator
    GeneRatioResult = Numerator / Denominator  # Perform the division
  )
# Sort the data frame by GeneRatioResult in descending order
top_go_results <- top_go_results %>%
  arrange(desc(GeneRatioResult))

# Create the plot using the reordered data
ggplot(top_go_results, aes(x = GeneRatioResult, y = reorder(Description, GeneRatioResult), size = Numerator, color = p.adjust)) +
  geom_point(alpha = 0.9) +  # Semi-transparent points
  theme_minimal() +  # Minimal theme for a clean look
  labs(
    x = "GeneRatio",
    y = "",
    size = "Gene count",
    color = "p.adjusted"
  ) +
  scale_color_gradient(high = "blue", low = "red") +  # Customize the color gradient
  scale_size_continuous(
    range = c(5, 15),  # Adjust the range for size (increase overall size of the circles)
    breaks = c(min(top_go_results$Numerator), median(top_go_results$Numerator), max(top_go_results$Numerator))
  ) +  # Only 3 bubbles in legend
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    plot.title = element_text(size = 20, face = "bold"),  # Increase title font size
    axis.title.x = element_text(size = 15),  # Increase x-axis title font size
    axis.title.y = element_text(size = 15),  # Increase y-axis title font size
    axis.text = element_text(size = 10, face = "bold"),  # Increase axis text font size
    legend.title = element_text(size = 15),  # Increase legend title font size
    legend.text = element_text(size = 12)  # Increase legend text font size
  )

