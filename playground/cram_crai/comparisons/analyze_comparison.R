#!/usr/bin/env Rscript

# Read the comparison data
data <- read.table("LPA_regions_comparison.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Remove rows where PERCENT_DIFF is N/A
data_clean <- data[data$PERCENT_DIFF != "N/A", ]
data_clean$PERCENT_DIFF <- as.numeric(data_clean$PERCENT_DIFF)

# Summary statistics
cat("=== LPA REGIONS COMPARISON SUMMARY ===\n")
cat("Number of samples analyzed:", nrow(data), "\n")
cat("Number of samples with coverage > 0:", nrow(data_clean), "\n\n")

cat("REGION 1 COVERAGE:\n")
cat("  Mean:", mean(data$REGION1_COVERAGE), "\n")
cat("  Median:", median(data$REGION1_COVERAGE), "\n")
cat("  Min:", min(data$REGION1_COVERAGE), "\n")
cat("  Max:", max(data$REGION1_COVERAGE), "\n\n")

cat("REGION 2 COVERAGE:\n")
cat("  Mean:", mean(data$REGION2_COVERAGE), "\n")
cat("  Median:", median(data$REGION2_COVERAGE), "\n")
cat("  Min:", min(data$REGION2_COVERAGE), "\n")
cat("  Max:", max(data$REGION2_COVERAGE), "\n\n")

if(nrow(data_clean) > 0) {
    cat("PERCENTAGE DIFFERENCE (Region1 vs Region2):\n")
    cat("  Mean:", mean(data_clean$PERCENT_DIFF), "%\n")
    cat("  Median:", median(data_clean$PERCENT_DIFF), "%\n")
    cat("  Min:", min(data_clean$PERCENT_DIFF), "%\n")
    cat("  Max:", max(data_clean$PERCENT_DIFF), "%\n")
    cat("  Std Dev:", sd(data_clean$PERCENT_DIFF), "%\n\n")
    
    # Correlation
    cor_result <- cor(data$REGION1_COVERAGE, data$REGION2_COVERAGE)
    cat("CORRELATION between regions:", cor_result, "\n\n")
    
    # Count samples with >10% difference
    high_diff <- sum(abs(data_clean$PERCENT_DIFF) > 10, na.rm=TRUE)
    cat("Samples with >10% difference:", high_diff, "out of", nrow(data_clean), "\n")
    
    # Samples with identical coverage
    identical_cov <- sum(data$REGION1_COVERAGE == data$REGION2_COVERAGE)
    cat("Samples with identical coverage:", identical_cov, "out of", nrow(data), "\n")
}

# Save detailed results
write.table(data, "LPA_comparison_detailed.txt", sep="\t", quote=FALSE, row.names=FALSE)
cat("\nDetailed results saved to: LPA_comparison_detailed.txt\n")
