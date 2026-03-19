# =============================================================================
# 02_query_neufer.R
# Purpose : Query all 37 MoTrPAC DA tables for the 85-gene Neufer framework
#           set. Join feature IDs → gene symbols via FEATURE_TO_GENE, filter
#           to the Neufer gene set, and export a unified long-format table.
# Inputs  : scripts/neufer_geneset.R   (sources the 85-gene set definition)
#           MotrpacRatTraining6moData   (all DA tables, FEATURE_TO_GENE)
# Outputs : data/neufer_da_all.csv      (42,770 rows × all assays/tissues)
#           data/neufer_da_significant.csv
#           data/nnt_da_all.csv
# Notes   : This script must be run before 03–07 — all downstream scripts
#           load data/neufer_da_all.csv rather than re-querying the package.
# =============================================================================
# 02: Query all DA tables for Neufer framework genes
# Strategy: join feature_ID from DA tables with FEATURE_TO_GENE
#   - TRNSCRPT: feature_ID = ENSRNOG... (Ensembl gene IDs)
#   - PROT/PHOSPHO/ACETYL/UBIQ: feature_ID = NP_/AP_... (protein accessions)
#   - METAB/IMMUNO: not gene-level, skipped
# ============================================================

library(MotrpacRatTraining6moData)
library(tidyverse)

source("scripts/neufer_geneset.R")

dir.create("data", showWarnings = FALSE)

# ----------------------------------------------------------
# 1. Build feature_ID → gene_symbol + framework category lookup
# ----------------------------------------------------------
cat("Building annotation lookup...\n")
data(FEATURE_TO_GENE)

# Create gene → category mapping from neufer_geneset
gene_category <- imap_dfr(neufer_geneset, function(genes, category) {
  tibble(gene_symbol = genes, framework_category = category)
}) %>%
  distinct(gene_symbol, .keep_all = TRUE)  # keep first category if gene appears in multiple

# Join annotation with Neufer gene set
neufer_annotation <- FEATURE_TO_GENE %>%
  filter(!is.na(gene_symbol), gene_symbol %in% all_neufer_genes) %>%
  select(feature_ID, gene_symbol, ensembl_gene, entrez_gene) %>%
  distinct() %>%
  left_join(gene_category, by = "gene_symbol")

cat("Neufer framework features in annotation:", nrow(neufer_annotation), "\n")
cat("Unique genes mapped:", n_distinct(neufer_annotation$gene_symbol), "of", length(all_neufer_genes), "\n\n")

# Genes not found in annotation
missing_genes <- setdiff(all_neufer_genes, neufer_annotation$gene_symbol)
if (length(missing_genes) > 0) {
  cat("Genes not found in FEATURE_TO_GENE:", paste(missing_genes, collapse = ", "), "\n\n")
}

# ----------------------------------------------------------
# 2. Discover all DA tables (skip METAB and IMMUNO)
# ----------------------------------------------------------
all_items <- data(package = "MotrpacRatTraining6moData")$results[, "Item"]
da_tables <- all_items[grepl("_DA$", all_items)]
da_tables <- da_tables[!grepl("^METAB_|^IMMUNO_|^ATAC_|^METHYL_", da_tables)]

cat("DA tables to query:", length(da_tables), "\n")
cat(paste(" ", da_tables), sep = "\n")
cat("\n")

# ----------------------------------------------------------
# 3. Query each DA table for Neufer genes
# ----------------------------------------------------------
query_da_table <- function(table_name) {
  # Load the table
  e <- new.env()
  data(list = table_name, package = "MotrpacRatTraining6moData", envir = e)
  da <- get(table_name, envir = e)

  # Join with Neufer annotation on feature_ID
  hits <- da %>%
    inner_join(neufer_annotation, by = "feature_ID") %>%
    select(
      gene_symbol, framework_category, feature_ID, ensembl_gene,
      assay, tissue, sex, comparison_group,
      logFC, adj_p_value, p_value,
      any_of(c("shrunk_logFC", "zscore", "selection_fdr", "tscore"))
    )

  if (nrow(hits) > 0) {
    cat(sprintf("  %-30s %d rows\n", table_name, nrow(hits)))
  }
  hits
}

cat("Querying DA tables...\n")
results_list <- map(da_tables, safely(query_da_table))

# Extract successful results
results <- map(results_list, "result") %>%
  compact() %>%
  bind_rows()

# Report errors
errors <- map2(da_tables, results_list, function(name, res) {
  if (!is.null(res$error)) cat("  ERROR in", name, ":", conditionMessage(res$error), "\n")
})

cat(sprintf("\nTotal results: %d rows across %d genes\n",
            nrow(results), n_distinct(results$gene_symbol)))

# ----------------------------------------------------------
# 4. NNT proof-of-concept summary
# ----------------------------------------------------------
cat("\n=== NNT across all tissues and assays ===\n")
nnt_results <- results %>%
  filter(gene_symbol == "Nnt") %>%
  arrange(adj_p_value, p_value)

cat("NNT rows:", nrow(nnt_results), "\n")
cat("Tissues with NNT data:", paste(unique(nnt_results$tissue), collapse = ", "), "\n")
cat("Assays with NNT data:", paste(unique(nnt_results$assay), collapse = ", "), "\n\n")

# NNT significant hits (adj_p < 0.05)
nnt_sig <- nnt_results %>% filter(adj_p_value < 0.05)
cat("NNT significant results (adj_p < 0.05):", nrow(nnt_sig), "\n")
if (nrow(nnt_sig) > 0) print(nnt_sig)

# NNT top results by p_value
cat("\nTop NNT results by p_value:\n")
print(nnt_results %>% select(assay, tissue, sex, comparison_group, logFC, adj_p_value, p_value) %>% head(20))

# ----------------------------------------------------------
# 5. Full gene set summary
# ----------------------------------------------------------
cat("\n=== Neufer gene set — significant hits (adj_p < 0.05) ===\n")
sig_hits <- results %>%
  filter(adj_p_value < 0.05) %>%
  arrange(framework_category, gene_symbol, adj_p_value)

cat("Total significant rows:", nrow(sig_hits), "\n")
cat("Significant genes:", n_distinct(sig_hits$gene_symbol), "\n\n")

sig_summary <- sig_hits %>%
  group_by(framework_category, gene_symbol, assay, tissue, sex) %>%
  summarise(
    n_timepoints_sig = n(),
    min_adj_p = min(adj_p_value),
    max_abs_logFC = max(abs(logFC)),
    .groups = "drop"
  ) %>%
  arrange(framework_category, min_adj_p)

print(sig_summary, n = 50)

# ----------------------------------------------------------
# 6. Save outputs
# ----------------------------------------------------------
cat("\nSaving results to data/...\n")

write_csv(results, "data/neufer_da_all.csv")
cat("  data/neufer_da_all.csv:", nrow(results), "rows\n")

write_csv(nnt_results, "data/nnt_da_all.csv")
cat("  data/nnt_da_all.csv:", nrow(nnt_results), "rows\n")

write_csv(sig_hits, "data/neufer_da_significant.csv")
cat("  data/neufer_da_significant.csv:", nrow(sig_hits), "rows\n")

write_csv(sig_summary, "data/neufer_da_significant_summary.csv")
cat("  data/neufer_da_significant_summary.csv:", nrow(sig_summary), "rows\n")

cat("\nDone.\n")
