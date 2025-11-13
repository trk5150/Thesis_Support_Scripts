#install.packages("openxlsx")
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(stringr)
library(ggplot2)
library(openxlsx)

# Load your data from Excel — assumes first column has UniProt IDs
data <- read.xlsx("C://Users//tik105//Downloads//Human Islet Pull-Down Data//4689 Human Islets Pull Down.xlsx")
uniprot_ids <- unique(data[[1]])  # Assumes first column contains UniProt IDs

# Convert UniProt IDs to Entrez IDs (clusterProfiler prefers Entrez)
converted_ids <- bitr(uniprot_ids,
                      fromType = "UNIPROT",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)

# Drop any NA conversions
entrez_ids <- unique(na.omit(converted_ids$ENTREZID))

# Define the GO enrichment function
perform_go_analysis <- function(genes, ont) {
  go_results <- enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ont,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  go_results_df <- as.data.frame(go_results)
  
  if (nrow(go_results_df) == 0) {
    message(paste("No significant GO terms found for", ont))
    return(NULL)
  }
  
  # Extract top terms and format for plotting
  top_go_results <- go_results_df %>%
    filter(Count > 3) %>%
    arrange(qvalue) %>%
    slice_head(n = 10) %>%     ##number of terms included on plots
    mutate(
      Description = ifelse(
        nchar(Description) > 70, 
        str_sub(Description, 1, str_locate(Description, ",")[, 1] - 1), 
        Description
      ),
      Numerator = as.numeric(str_extract(GeneRatio, "^[^/]+")),
      Denominator = as.numeric(str_extract(GeneRatio, "(?<=/)[0-9]+")),
      GeneRatioResult = Numerator / Denominator
    )
  
  # Bubble plot
  plot <- ggplot(top_go_results, aes(x = qvalue, y = reorder(Description, -qvalue), size = Numerator, color = GeneRatioResult)) +
    geom_point(alpha = 0.9) +
    theme_minimal() +
    labs(
      title = paste("GO Analysis -", ont),
      x = "q-value",
      y = "",
      size = "Gene count",
      color = "Gene Ratio"
    ) +
    scale_x_reverse() +
    scale_color_gradient(low = "blue", high = "red") +
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
  
  # Save plot and data
  output_directory <- "C://Users//tik105//Desktop//Pulldown go enrichment//"
  dir.create(output_directory, showWarnings = FALSE)
  
  plot_file <- paste0(output_directory, "GO_Plot_", ont, ".png")
  ggsave(plot_file, plot, width = 10, height = 6)
  
  excel_file <- paste0(output_directory, "GO_Results_", ont, ".xlsx")
  write.xlsx(go_results_df, file = excel_file)
}

# Loop through GO categories
for (ont in c("BP", "MF", "CC")) { ##Can contain "BP" "MF" or "CC" and any combination of those, separated by ","
  perform_go_analysis(entrez_ids, ont)
}

