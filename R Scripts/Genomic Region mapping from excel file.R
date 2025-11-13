library(ggplot2)
library(readxl)
library(dplyr)

# Load data
file_path <- "C://Users//tik105//OneDrive - Harvard University//Atollin Paper//Genomic Diagraming//Mouse coordinates.xlsx"  # Change to your file path
#file_path <- "C://Users//tik105//OneDrive - Harvard University//Atollin Paper//Genomic Diagraming//Human coordinates.xlsx"  # Change to your file path
df <- read_excel(file_path)

# Reduce heights for all rectangles
exon_height <- 0.01
non_coding_height <- 0.005
feature_height <- 0.01

# Remove overlap between "Exon" and "Non-coding exon"
df_exon <- df %>% filter(Type == "Exon")
df_nc_exon <- df %>% filter(Type == "Non-coding exon")

# Adjust exons by removing overlapping regions with "Non-coding exon"
for (i in seq_len(nrow(df_nc_exon))) {
  nc_start <- df_nc_exon$Start[i]
  nc_end <- df_nc_exon$End[i]
  nc_name <- df_nc_exon$Name[i]
  
  df_exon <- df_exon %>%
    rowwise() %>%
    mutate(
      Start = ifelse(Name == nc_name & Start < nc_end & End > nc_start, max(Start, nc_end), Start),
      End = ifelse(Name == nc_name & Start < nc_end & End > nc_start, min(End, nc_start), End)
    ) %>%
    ungroup()
}

# Filter out exons that got fully removed by overlap
df_exon <- df_exon %>% filter(Start < End)

# Combine exon and non-coding exon for label preservation
df_labels <- df_exon %>%
  bind_rows(df_nc_exon) %>%
  distinct(Name, .keep_all = TRUE)  # Ensure no duplicate labels

# Offset overlapping feature elements
df_features <- df %>% filter(Type == "Feature") %>%
  arrange(Start) %>%
  mutate(feature_offset = seq(0.1, 0.2, length.out = n()))  # Stagger features

# Create the plot
ggplot() +
  
  # Plot "Total" as a black line
  geom_segment(data = df %>% filter(Type == "Total"),
               aes(x = Start, xend = End, y = 0, yend = 0),
               color = "black", linewidth = 1) +  
  
  # Plot adjusted "Exon" as purple rectangles
  geom_rect(data = df_exon,
            aes(xmin = Start, xmax = End, ymin = -exon_height, ymax = exon_height, fill = "Exon"),
            color = "black") +
  
  # Plot "Non-coding exon" as smaller purple rectangles
  geom_rect(data = df_nc_exon,
            aes(xmin = Start, xmax = End, ymin = -non_coding_height, ymax = non_coding_height, fill = "Non-coding exon"),
            color = "black") +
  
  # Write exon labels **vertically**
  geom_text(data = df_labels,
            aes(x = (Start + End) / 2, y = -exon_height - 0.02, label = Name), 
            size = 3, angle = 90) +  # **Rotated text**
  
  # Plot "Feature" as orange rectangles with staggered heights
  geom_rect(data = df_features,
            aes(xmin = Start, xmax = End, ymin = feature_offset, ymax = feature_offset + feature_height, fill = "Feature"),
            color = "black") +
  
  # Add staggered labels for Features (horizontal)
  geom_text(data = df_features,
            aes(x = (Start + End) / 2, y = feature_offset + feature_height + 0.005, label = Name), 
            size = 3, vjust = 0) +
  
  # Define colors for fill
  scale_fill_manual(values = c("Exon" = "purple", "Non-coding exon" = "mediumpurple", "Feature" = "black")) +  
  
  # Customize the theme
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  )+
  
  labs(title = "Mouse Genomic Features Diagram", x = "Genomic Coordinates")
