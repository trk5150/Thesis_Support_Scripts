#!/bin/R
# Protein Plotting Script for Two Architecture Files with "Miss" Exclusion, Custom Black Bar, and Updated File Paths

# Install required packages if not already installed
#if (!require("seqinr")) install.packages("seqinr")

library("seqinr")

# --- User-defined parameters ---
proteinArchitectureFile1 <- "C://Users//tik105//OneDrive - Harvard University//Atollin Paper//Conservation//PlotProtein//C1orf127_architecture_file.txt"  # Path to the first architecture file
proteinArchitectureFile2 <- "C://Users//tik105//OneDrive - Harvard University//Atollin Paper//Conservation//PlotProtein//C1orf127_architecture_file - Short isoform.txt"  # Path to the second architecture file
alignmentFile <- "C://Users//tik105//OneDrive - Harvard University//Atollin Paper//Conservation//Aligned_Model-Orgs_sequence.fas"
referenceSequencePositionInFile <- 1                # Reference sequence position in the alignment
nameOfYourQuery <- "C1orf127 subset"                    # Name of the protein for the plot
tickSize <- 10                                      # Tick interval for the x-axis

# --- Load input files ---
# Read the architecture files
pa1 <- read.table(proteinArchitectureFile1, sep="\t", header=TRUE)
pa2 <- read.table(proteinArchitectureFile2, sep="\t", header=TRUE)

# Identify the domain labeled "Miss" in the second architecture file
miss_domain <- subset(pa2, architecture_name == "Miss")
pa2 <- subset(pa2, architecture_name != "Miss")  # Exclude "Miss" from plotting

# Determine the starting position for the black bar
black_bar_start <- if (nrow(miss_domain) > 0) {
  max(miss_domain$end_site)  # End of the "Miss" domain
} else {
  0  # Default to 0 if no "Miss" domain is present
}

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
tiff(paste(nameOfYourQuery, "_protein_plot.tif", sep=""), height=8.5, width=11, units="in", res=300)
layout(matrix(c(1, 2), nrow=1), widths=c(1, 3))
par(oma=c(4, 0, 4, 0), mar=c(5, 0, 4, 0) + 0.4)

# Legend
plot((-30:-15), rep(-1, 16), col="purple3", type="l", ann=FALSE, bty="n", xaxt="n", yaxt="n", xlim=c(-160, -15), ylim=c(1, -8))
lines((-30:-15), rep(0, 16), col="purple3")
lines((-30:-15), rep(-0.5, 16), col="purple3")
text(-60, -0.5, "Conservation", col="purple3", cex=0.9, font=2, srt=90)
text(-45, -1, "1", col="purple3", cex=0.9)
text(-45, -0.5, "0.5", col="purple3", cex=0.9)
text(-45, 0, "0", col="purple3", cex=0.9)

# Main plot
plot((1:length(scorePlot)), rep(-2, length(scorePlot)), type="l", lwd=5, main=paste("Protein Architecture and Conservation of", nameOfYourQuery),
     xlab="Amino Acid Position", ylab="", xlim=c(0, length(scorePlot)), ylim=c(1, -8), cex.lab=0.9, cex.main=1, yaxt="n", bty="n", font=2, xaxt="n")

# Conservation score line
lines(scorePlot, col="purple3")

# Ticks on x-axis
ticks <- seq(0, length(scorePlot), by=tickSize)
axis(side=1, at=ticks, las=3)

# Plot architecture from the first file
for (i in 1:nrow(pa1)) {
  rect(pa1$start_site[i], -2.1, pa1$end_site[i], -1.9, col="lightseagreen", border="black")
  text(median(c(pa1$start_site[i], pa1$end_site[i])), -1.7, pa1$architecture_name[i], cex=1, col="darkblue")
}

# Add a centered black bar for the second architecture plot
black_bar_center <- mean(c(-4.3, -4.1))  # Center vertically
segments(black_bar_start, black_bar_center, length(scorePlot), black_bar_center, col="black", lwd=5)

# Plot architecture from the second file (excluding "Miss")
for (i in 1:nrow(pa2)) {
  rect(pa2$start_site[i], -4.3, pa2$end_site[i], -4.1, col="lightcoral", border="black")
  text(median(c(pa2$start_site[i], pa2$end_site[i])), -3.9, pa2$architecture_name[i], cex=1, col="darkred")
}

dev.off()

# Inform the user of output
cat("Plot saved as:", paste(nameOfYourQuery, "_protein_plot.pdf", sep=""), "\n")
