# =============================================================================
# 01_explore_data.R
# Purpose : Explore the MoTrPAC data structure — available DA tables, assay
#           types, tissue codes, feature ID formats, and gene annotation.
# Inputs  : MotrpacRatTraining6moData R package (auto-loaded)
# Outputs : console output only (no saved files)
# Notes   : Run this first to understand naming conventions used by all
#           downstream scripts. No files are written.
# =============================================================================

library(MotrpacRatTraining6mo)
library(MotrpacRatTraining6moData)
library(tidyverse)

# What data objects exist?
available <- data(package = "MotrpacRatTraining6moData")
all_items <- available$results[, "Item"]
cat("Total data objects:", length(all_items), "\n\n")

# Separate DA tables from normalized data
da_items <- all_items[grepl("_DA$", all_items)]
norm_items <- all_items[grepl("NORM_DATA$|_NORM$", all_items)]

cat("Differential analysis tables:\n")
print(da_items)
cat("\nNormalized data tables:\n")
print(norm_items)

# Load skeletal muscle transcriptomics as our reference
data(TRNSCRPT_SKMGN_DA)
cat("\n\n=== TRNSCRPT_SKMGN_DA structure ===\n")
cat("Dimensions:", dim(TRNSCRPT_SKMGN_DA), "\n")
cat("Columns:", paste(colnames(TRNSCRPT_SKMGN_DA), collapse = ", "), "\n")
cat("\nFirst rows:\n")
print(head(TRNSCRPT_SKMGN_DA, 10))

cat("\nComparison groups:\n")
if ("comparison_group" %in% colnames(TRNSCRPT_SKMGN_DA)) {
  print(unique(TRNSCRPT_SKMGN_DA$comparison_group))
}

# Quick check: is NNT in the data?
cat("\n=== NNT check ===\n")
nnt_check <- TRNSCRPT_SKMGN_DA[grepl("Nnt|NNT", TRNSCRPT_SKMGN_DA$feature_ID, ignore.case = TRUE), ]
if (nrow(nnt_check) > 0) {
  print(nnt_check)
} else {
  cat("NNT not found by feature_ID, checking other columns...\n")
  for (col in colnames(TRNSCRPT_SKMGN_DA)) {
    if (is.character(TRNSCRPT_SKMGN_DA[[col]])) {
      hits <- TRNSCRPT_SKMGN_DA[grepl("Nnt", TRNSCRPT_SKMGN_DA[[col]], ignore.case = TRUE), ]
      if (nrow(hits) > 0) {
        cat("Found NNT in column:", col, "\n")
        print(hits)
      }
    }
  }
}

cat("\n=== Setup complete. Ready for analysis. ===\n")
