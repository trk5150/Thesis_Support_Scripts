##Code authored by Tim Kunz,
##Starting from a raw counts matrix from bulk RNA-seq data, code performs differential expression analysis
##then performs GO:BP and GO:MF analysis for upregulated and downregulated genes. Saves excel summary of GO analyses and bubble plots
##thresholds, labels and file locations are not taken as arguments; script contents must be editted



##Code for installing necessary packages; there may be additional required packages 
# Install BiocManager if it's not installed
#if (!requireNamespace("BiocManager", quietly = TRUE)) {
#  install.packages("BiocManager")
#}

# Install necessary Bioconductor packages
#BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "DOSE"), ask = FALSE)

# Load required libraries
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
library(writexl)

# Load your data
data <- read_excel("C://Users//tik105//Desktop//Neonate phenotyping 5-1-2024//Bulk RNA-seq//Null pancreas//Results//Counts.xlsx")


# Prepare count data
#reads in the count matrix from the above excel file, specifically from the following columns
count_data <- as.matrix(data[, c("7-1", "7-2", "7-5", "13-1", "13-2", "13-5", "13-7")]) #these are animal identifiers (litter#-pup#)
rownames(count_data) <- data$Gene_Name

# Prepare metadata
#assigns experimental condition to each above animal
coldata <- data.frame(
  condition = c("Het", "Het", "Null", "Het", "Null", "Het", "Null"),
  row.names = colnames(count_data)
)
coldata$condition <- factor(coldata$condition)

# DESeq2 Analysis
#DESeq2 is the standrard RNA-seq analysis tools for bulk sequencing data
#performs normalization, statistics and corrections
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = coldata, design = ~ condition)
dds <- DESeq(dds)


# Get differential expression results
res <- results(dds, contrast = c("condition", "Null", "Het"))
res <- na.omit(res)  # Remove rows with NA values
res <- res[order(res$padj), ]  # Sort by adjusted p-value

#From DESeq2 results, extraction information for GO:analysis
#Extract fold changes
fold_changes <- res$log2FoldChange
names(fold_changes) <- rownames(res)

#Set thresholds for genes considered for GO:enrichment
#can tighten thresholds if too many genes are being included, but shouldn't loosen them much
# Define thresholds
alpha <- 0.001  # Adjusted p-value threshold (very low p-value)
lfc_threshold <- 0.5  # Log2 fold change threshold (~2x average expression change)

# Filter upregulated and downregulated genes
upregulated_genes <- rownames(res[res$padj < alpha & res$log2FoldChange > lfc_threshold, ])
downregulated_genes <- rownames(res[res$padj < alpha & res$log2FoldChange < -lfc_threshold, ])

#Defines the GO:enrichment function
# Function to perform GO enrichment and save results
perform_go_analysis <- function(genes, direction, ont) {
  go_results <- enrichGO(
    gene = genes,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = ont,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  # Convert results to data frame
  go_results_df <- as.data.frame(go_results)
  
  # Extract top GO terms
  top_go_results <- go_results_df %>%
    filter(Count > 3) %>%                   # Only include terms with more than 3 genes
    arrange(qvalue) %>%                     # Sort by qvalue (lowest first)
    slice_head(n = 10) %>%                  # Select the top 10 terms
    mutate(
      Description = ifelse(
        nchar(Description) > 70, 
        str_sub(Description, 1, str_locate(Description, ",")[, 1] - 1), 
        Description
      ),
      Numerator = as.numeric(str_extract(GeneRatio, "^[^/]+")),  # Extract numerator
      Denominator = as.numeric(str_extract(GeneRatio, "(?<=/)[0-9]+")),  # Extract denominator
      GeneRatioResult = Numerator / Denominator  # Perform the division
    )
  #Generate a bubble plot for each analysis
  # Plot
  plot <- ggplot(top_go_results, aes(x = qvalue, y = reorder(Description, -qvalue), size = Numerator, color = GeneRatioResult)) +
    geom_point(alpha = 0.9) +
    theme_minimal() +
    labs(
      title = paste("GO Analysis -", direction, "-", ont),
      x = "q-value",
      y = "",
      size = "Gene count",
      color = "Gene Ratio"
    ) +
    scale_x_reverse() +  # Reverse the x-axis
    scale_color_gradient(low = "blue", high = "red") +  # Customize the color gradient
    scale_size_continuous(
      range = c(5, 9),
      breaks = c(min(top_go_results$Numerator), max(top_go_results$Numerator))
    ) +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      plot.title = element_text(size = 20, face = "bold"),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      axis.text = element_text(size = 10, face = "bold"),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 12)
    )
  
  #Output directory where plots and excel files will be saved
  #If the plots automatically generated don't look exactly as you'd like, can use the saved excel files to generate custom plots
  output_directory <- "C://Users//tik105//Desktop//Neonate phenotyping 5-1-2024//Bulk RNA-seq//Null pancreas//Results//"
  # Save plot
  plot_file <- paste0(output_directory,"GO_Plot_", direction, "_", ont, ".png")
  print(plot_file)
  ggsave(plot_file, plot)
  
  # Save results
  excel_file <- paste0(output_directory,"GO_Results_", direction, "_", ont, ".xlsx")
  print(excel_file)
  write_xlsx(go_results_df, path = excel_file)
}

#Loop to perform go analysis for pairwise combinations of biological process and molecular function for genes upregulated and genes downregulated
# Perform GO analysis for combinations
for (direction in c("Upregulated", "Downregulated")) {
  for (ont in c("BP", "MF")) {
    genes <- if (direction == "Upregulated") upregulated_genes else downregulated_genes
    perform_go_analysis(genes, direction, ont)
  }
}
