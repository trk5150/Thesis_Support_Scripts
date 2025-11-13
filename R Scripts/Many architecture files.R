#!/bin/R
# Protein Plotting Script for Four Architecture Files with "Miss" Exclusion and Conservation Trace

library("seqinr")

# --- User-defined parameters ---
proteinArchitectureFile3 <- "C://Users//tik105//OneDrive - Harvard University//Atollin Paper//Conservation//PlotProtein//C1orf127_architecture_file.txt"
proteinArchitectureFile4 <- "C://Users//tik105//OneDrive - Harvard University//Atollin Paper//Conservation//PlotProtein//C1orf127_architecture_file - Short isoform.txt"
#proteinArchitectureFile2 <- "C://Users//tik105//OneDrive - Harvard University//Atollin Paper//Conservation//PlotProtein//C1orf127_architecture_file - Z domains.txt"
proteinArchitectureFile1 <- "C://Users//tik105//OneDrive - Harvard University//Atollin Paper//Conservation//PlotProtein//C1orf127_architecture_file - Exons.txt"
proteinArchitectureFile2 <- "C://Users//tik105//OneDrive - Harvard University//Atollin Paper//Conservation//PlotProtein//C1orf127_architecture_file - 13.5 del.txt"

alignmentFile <- "C://Users//tik105//OneDrive - Harvard University//Atollin Paper//Conservation//Aligned_Model-Orgs_sequence.fas"

referenceSequencePositionInFile <- 1
nameOfYourQuery <- "C1orf127 subset"
tickSize <- 10

# --- Load input files ---
pa1 <- read.table(proteinArchitectureFile1, sep="\t", header=TRUE)
pa2 <- read.table(proteinArchitectureFile2, sep="\t", header=TRUE)
pa3 <- read.table(proteinArchitectureFile3, sep="\t", header=TRUE)
pa4 <- read.table(proteinArchitectureFile4, sep="\t", header=TRUE)

# Identify "Miss" regions and remove them from plotting
miss1 <- subset(pa1, architecture_name == "Miss")
miss2 <- subset(pa2, architecture_name == "Miss")
miss3 <- subset(pa3, architecture_name == "Miss")
miss4 <- subset(pa4, architecture_name == "Miss")

pa1 <- subset(pa1, architecture_name != "Miss")
pa2 <- subset(pa2, architecture_name != "Miss")
pa3 <- subset(pa3, architecture_name != "Miss")
pa4 <- subset(pa4, architecture_name != "Miss")

# Read the alignment file
a <- read.fasta(alignmentFile)

# --- Prepare conservation data ---
seq <- list()
for (i in 1:length(a)) {
  seq[[i]] <- a[[i]][1:length(a[[i]])]
}

numberOfSeq <- length(seq)
mat <- matrix(0, nrow=length(a), ncol=length(a[[1]]))
for (i in 1:length(seq)) {
  mat[i, ] <- seq[[i]]
}

df <- as.data.frame(mat)
tdf <- t(df)
referenceSeq <- tdf[which(tdf[, as.numeric(referenceSequencePositionInFile)] != "-"), ]
referenceSeq <- as.data.frame(referenceSeq)
write.table(referenceSeq, file="alignment_table", sep="\t", quote=F, row.names=F, col.names=F)

# Calculate conservation scores
a <- read.table("alignment_table", sep="\t")
a <- data.frame(lapply(a, as.character), stringsAsFactors=FALSE)

for (i in 1:nrow(a)) {
  a[i, "consensus"] <- paste(as.character(a[i, ]), collapse="")
}

countBases <- function(string) {
  table(strsplit(string, "")[[1]])
}

c <- as.character(a[, "consensus"])
tab <- list()
for (i in 1:length(c)) {
  tab[[i]] <- countBases(c[i])
}

score <- rep(0, nrow(a))
for (i in 1:length(tab)) {
  for (j in 1:length(tab[[i]])) {
    if ((names(tab[[i]][j])) == a[i, ][as.numeric(referenceSequencePositionInFile)])
      score[i] <- tab[[i]][j]
  }
}
scorePlot <- -(((score / numberOfSeq)))
system("rm alignment_table")

# --- Plotting ---
tiff(paste(nameOfYourQuery, "_protein_plot.tif", sep=""), height=10, width=11, units="in", res=300)
layout(matrix(c(1, 2), nrow=1), widths=c(1, 3))
par(oma=c(4, 0, 4, 0), mar=c(5, 0, 4, 0) + 0.4)

# Legend
plot((-30:-15), rep(-1, 16), col="purple3", type="l", ann=FALSE, bty="n", xaxt="n", yaxt="n", xlim=c(-160, -15), ylim=c(1, -10))
lines((-30:-15), rep(0, 16), col="purple3")
lines((-30:-15), rep(-0.5, 16), col="purple3")
text(-60, -0.5, "Conservation", col="purple3", cex=0.9, font=2, srt=90)
text(-45, -1, "1", col="purple3", cex=0.9)
text(-45, -0.5, "0.5", col="purple3", cex=0.9)
text(-45, 0, "0", col="purple3", cex=0.9)

# Main plot
plot((1:length(scorePlot)), type="n", xlab="Amino Acid Position", ylab="Conservation Score",
     xlim=c(0, length(scorePlot)), ylim=c(1, -10), cex.lab=0.9, cex.main=1, yaxt="n", xaxt="n", bty="n",
     main=paste("Protein Architecture and Conservation of", nameOfYourQuery))

# Conservation score line
lines(scorePlot, col="purple3")

# Ticks on x-axis
ticks <- seq(0, length(scorePlot), by=tickSize)
axis(side=1, at=ticks, las=3, line=-3.3)  # Moves x-axis higher

# Define architecture levels and black bar positions
architecture_levels <- c(-2, -4, -6, -8)
black_bar_levels <- architecture_levels - 0.2  

# Function to plot architecture and black bars while omitting "Miss" regions
plot_architecture_with_black_bar <- function(pa, miss, level, color, text_color, black_bar_level) {
  # Draw black bar excluding "Miss" regions
  if (nrow(miss) > 0) {
    for (i in 1:nrow(miss)) {
      if (i == 1 && miss$start_site[i] > 1) {
        segments(0, black_bar_level, miss$start_site[i] - 1, black_bar_level, col="black", lwd=5)
      }
      if (i < nrow(miss)) {
        segments(miss$end_site[i], black_bar_level, miss$start_site[i+1] - 1, black_bar_level, col="black", lwd=5)
      } else {
        segments(miss$end_site[i], black_bar_level, length(scorePlot), black_bar_level, col="black", lwd=5)
      }
    }
  } else {
    segments(0, black_bar_level, length(scorePlot), black_bar_level, col="black", lwd=5)
  }
  
  # Draw architecture elements
  for (i in 1:nrow(pa)) {
    rect(pa$start_site[i], level - 0.3, pa$end_site[i], level - 0.1, col=color, border="black")
    text(median(c(pa$start_site[i], pa$end_site[i])), level + 0.1, pa$architecture_name[i], cex=1, col=text_color)
  }
}

# Plot architectures and conservation trace
plot_architecture_with_black_bar(pa1, miss1, architecture_levels[1], "lightseagreen", "darkblue", black_bar_levels[1])
plot_architecture_with_black_bar(pa2, miss2, architecture_levels[2], "black", "darkblue", black_bar_levels[2])
plot_architecture_with_black_bar(pa3, miss3, architecture_levels[3], "lightblue", "darkblue", black_bar_levels[3])
plot_architecture_with_black_bar(pa4, miss4, architecture_levels[4], "purple", "darkblue", black_bar_levels[4])

dev.off()

# Inform the user of output
cat("Plot saved as:", paste(nameOfYourQuery, "_protein_plot.tif", sep=""), "\n")
