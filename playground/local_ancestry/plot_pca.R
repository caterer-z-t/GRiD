#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
cat("📌 Command-line arguments:\n")
print(args)

if (length(args) < 2) {
  stop("❌ You must provide two PC numbers (e.g., 1 2 for PC1 vs PC2).")
}

PC_A <- as.integer(args[1])
PC_B <- as.integer(args[2])

# Define input file prefix
prefix <- "merged_SOL_ref"

# Load eigenvectors
evec <- read_table2(paste0(prefix, ".eigenvec"), col_names = FALSE)
colnames(evec)[1:2] <- c("FID", "IID")
colnames(evec)[3:ncol(evec)] <- paste0("PC", 1:(ncol(evec) - 2))

# Load eigenvalues and compute variance explained
eval <- scan(paste0(prefix, ".eigenval"))
var_explained <- round(100 * eval / sum(eval), 2)

# Check if requested PCs are valid
num_pcs <- length(eval)
if (PC_A > num_pcs || PC_B > num_pcs) {
  stop(paste0("❌ Requested PC index exceeds number of PCs available (", num_pcs, ")."))
}

# Dynamically build column names
x_col <- paste0("PC", PC_A)
y_col <- paste0("PC", PC_B)

# Prepare plot data
plot_df <- evec %>%
  select(IID, all_of(x_col), all_of(y_col))

# Make plot
p <- ggplot(plot_df, aes_string(x = x_col, y = y_col, label = "IID")) +
  geom_point(alpha = 0.7, color = "#377eb8", size = 1) +
  # geom_text(size = 2.5, vjust = -0.5, check_overlap = TRUE) +
  theme_minimal(base_size = 14) +
  labs(
    title = paste("PCA:", x_col, "vs", y_col),
    x = paste0(x_col, " (", var_explained[PC_A], "% variance)"),
    y = paste0(y_col, " (", var_explained[PC_B], "% variance)")
  )

# Save plot
output_file <- paste0(prefix, "_", x_col, "_", y_col, ".png")
ggsave(filename = paste0("figures/", output_file), plot = p, width = 8, height = 6, dpi = 300)

cat("✅ PCA plot saved as", output_file, "\n")
