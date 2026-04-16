# =============================================================================
# explore_aging_signals.R
# Purpose : Reconnaissance — check whether canonical aging biomarker genes
#           show any signal in MoTrPAC rat endurance training data.
#           Also inventory METHYL DA tables to assess epigenetic clock feasibility.
# Inputs  : MotrpacRatTraining6moData (DA tables, FEATURE_TO_GENE)
# Outputs : data/aging_recon.csv (flat table for eyeballing)
#           console summary
# Separate from the Neufer redox framework — exploratory only.
# =============================================================================

library(MotrpacRatTraining6moData)
library(tidyverse)

dir.create("data", showWarnings = FALSE)

# ----------------------------------------------------------
# 1. Define ~15 canonical aging biomarker genes
#    Source: Wu et al. (2025) Nat Rev Mol Cell Biol 26:826-847
#    NOT a framework — just the most-cited markers from the review
# ----------------------------------------------------------
aging_genes <- list(
  # Cell senescence
  senescence = c("Cdkn2a",    # p16^INK4a
                 "Cdkn1a",    # p21^CIP1
                 "Trp53",     # p53
                 "Rb1"),      # Rb

  # SASP / inflammaging
  sasp = c("Gdf15",     # growth differentiation factor 15 — conserved primate marker
           "Serpine1",   # PAI-1 (plasminogen activator inhibitor)
           "Il6",        # interleukin-6
           "Ccl2",       # MCP-1
           "Tnf"),       # TNF-alpha

  # Nuclear architecture
  nuclear = c("Lmnb1",     # lamin B1 (lost with age)
              "Lmnb2",     # lamin B2
              "Lmna"),     # lamin A/C

  # NAD+ metabolism
  nad = c("Nampt",     # rate-limiting NAD+ salvage
          "Nmnat1",    # NAD+ synthesis
          "Nmnat2",
          "Sirt1",     # NAD+-dependent deacetylase (nuclear)
          "Sirt3",     # NAD+-dependent deacetylase (mitochondrial)
          "Sirt6"),    # NAD+-dependent deacetylase (chromatin)

  # Epigenetic machinery (writers/erasers of aging-associated marks)
  epigenetic = c("Dnmt1",     # DNA methyltransferase 1 (maintenance)
                 "Dnmt3a",    # DNA methyltransferase 3a
                 "Tet2",      # ten-eleven translocation 2
                 "Suv39h1",   # H3K9me3 methyltransferase
                 "Ehmt2",     # G9a — H3K9me2
                 "Kdm4a"),    # H3K9 demethylase

  # Lysosomal / autophagy
  lysosomal = c("Glb1",      # beta-galactosidase (SA-beta-gal)
                "Lamp1",     # lysosome marker
                "Tfeb",      # master autophagy regulator
                "Map1lc3b"), # LC3B — autophagosome marker

  # Telomere-associated
  telomere = c("Tert",      # telomerase reverse transcriptase
               "Terf2"),    # shelterin component

  # Mitophagy / mito quality
  mitophagy = c("Pink1",     # PTEN-induced kinase 1
                "Prkn",      # parkin
                "Bnip3l",    # Nix — mitophagy receptor
                "Fundc1"),   # mitophagy receptor

  # cGAS-STING pathway (links DNA damage → inflammation)
  cgas_sting = c("Cgas",      # cyclic GAS
                 "Sting1",    # STING (Tmem173)
                 "Mb21d1")    # alternate name for cGAS in some species
)

all_aging_genes <- unique(unlist(aging_genes))
cat("Aging biomarker probe set:", length(all_aging_genes), "genes\n\n")

# ----------------------------------------------------------
# 2. Build feature_ID → gene_symbol lookup (same as 02_query)
# ----------------------------------------------------------
cat("Building annotation lookup...\n")
data(FEATURE_TO_GENE)

gene_category <- imap_dfr(aging_genes, function(genes, category) {
  tibble(gene_symbol = genes, aging_category = category)
}) %>%
  distinct(gene_symbol, .keep_all = TRUE)

aging_annotation <- FEATURE_TO_GENE %>%
  filter(!is.na(gene_symbol), gene_symbol %in% all_aging_genes) %>%
  select(feature_ID, gene_symbol, ensembl_gene, entrez_gene) %>%
  distinct() %>%
  left_join(gene_category, by = "gene_symbol")

cat("Aging genes in FEATURE_TO_GENE:", n_distinct(aging_annotation$gene_symbol),
    "of", length(all_aging_genes), "\n")

missing <- setdiff(all_aging_genes, aging_annotation$gene_symbol)
if (length(missing) > 0) {
  cat("NOT found in FEATURE_TO_GENE:", paste(missing, collapse = ", "), "\n")
}
cat("\n")

# ----------------------------------------------------------
# 3. Discover DA tables (include same assays as 02_query)
# ----------------------------------------------------------
all_items <- data(package = "MotrpacRatTraining6moData")$results[, "Item"]
da_tables <- all_items[grepl("_DA$", all_items)]
gene_da_tables <- da_tables[!grepl("^METAB_|^IMMUNO_|^ATAC_|^METHYL_", da_tables)]

cat("Gene-level DA tables:", length(gene_da_tables), "\n\n")

# ----------------------------------------------------------
# 4. Query each DA table
# ----------------------------------------------------------
query_da_table <- function(table_name) {
  e <- new.env()
  data(list = table_name, package = "MotrpacRatTraining6moData", envir = e)
  da <- get(table_name, envir = e)

  hits <- da %>%
    inner_join(aging_annotation, by = "feature_ID") %>%
    select(
      gene_symbol, aging_category, feature_ID, ensembl_gene,
      assay, tissue, sex, comparison_group,
      logFC, adj_p_value, p_value,
      any_of(c("shrunk_logFC", "zscore", "selection_fdr", "tscore"))
    )

  if (nrow(hits) > 0) {
    cat(sprintf("  %-30s %d rows\n", table_name, nrow(hits)))
  }
  hits
}

cat("Querying DA tables for aging genes...\n")
results_list <- map(gene_da_tables, safely(query_da_table))

results <- map(results_list, "result") %>%
  compact() %>%
  bind_rows()

errors <- map2(gene_da_tables, results_list, function(name, res) {
  if (!is.null(res$error)) cat("  ERROR in", name, ":", conditionMessage(res$error), "\n")
})

cat(sprintf("\nTotal: %d rows across %d genes\n\n",
            nrow(results), n_distinct(results$gene_symbol)))

# ----------------------------------------------------------
# 5. Summary: significant hits per gene at 8w
# ----------------------------------------------------------
cat("=== Significant hits (adj_p < 0.05, 8w training) ===\n")
sig_8w <- results %>%
  filter(adj_p_value < 0.05, comparison_group == "8w") %>%
  count(aging_category, gene_symbol, assay, tissue, sex, name = "n_features") %>%
  arrange(aging_category, gene_symbol)

if (nrow(sig_8w) > 0) {
  print(as.data.frame(sig_8w), row.names = FALSE)
} else {
  cat("(none)\n")
}

cat("\n=== Significant hits (adj_p < 0.05, any timepoint) ===\n")
sig_any <- results %>%
  filter(adj_p_value < 0.05) %>%
  count(aging_category, gene_symbol, name = "n_sig_rows") %>%
  arrange(desc(n_sig_rows))

if (nrow(sig_any) > 0) {
  print(as.data.frame(sig_any), row.names = FALSE)
} else {
  cat("(none)\n")
}

# ----------------------------------------------------------
# 6. Top hits by effect size at 8w
# ----------------------------------------------------------
cat("\n=== Largest fold changes at 8w (top 30) ===\n")
top_fc <- results %>%
  filter(comparison_group == "8w") %>%
  arrange(desc(abs(logFC))) %>%
  select(gene_symbol, aging_category, assay, tissue, sex, logFC, adj_p_value) %>%
  head(30)
print(as.data.frame(top_fc), row.names = FALSE)

# ----------------------------------------------------------
# 7. Per-gene summary across all conditions
# ----------------------------------------------------------
cat("\n=== Per-gene summary (all conditions) ===\n")
gene_summary <- results %>%
  group_by(aging_category, gene_symbol) %>%
  summarise(
    n_rows = n(),
    n_tissues = n_distinct(tissue),
    n_assays = n_distinct(assay),
    n_sig_005 = sum(adj_p_value < 0.05, na.rm = TRUE),
    n_sig_8w = sum(adj_p_value < 0.05 & comparison_group == "8w", na.rm = TRUE),
    max_abs_logFC = max(abs(logFC), na.rm = TRUE),
    mean_logFC_8w = mean(logFC[comparison_group == "8w"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(aging_category, desc(n_sig_005))

print(as.data.frame(gene_summary), row.names = FALSE)

# ----------------------------------------------------------
# 8. METHYL assay inventory
# ----------------------------------------------------------
cat("\n\n========================================\n")
cat("=== METHYL DA TABLE INVENTORY ===\n")
cat("========================================\n\n")

methyl_tables <- da_tables[grepl("^METHYL_", da_tables)]
cat("METHYL DA tables:", length(methyl_tables), "\n")
cat(paste(" ", methyl_tables), sep = "\n")
cat("\n")

if (length(methyl_tables) > 0) {
  for (tbl in methyl_tables) {
    e <- new.env()
    data(list = tbl, package = "MotrpacRatTraining6moData", envir = e)
    da <- get(tbl, envir = e)

    cat(sprintf("\n--- %s ---\n", tbl))
    cat("Dimensions:", nrow(da), "rows x", ncol(da), "cols\n")
    cat("Columns:", paste(names(da), collapse = ", "), "\n")
    cat("Sample feature_IDs:\n")
    print(head(unique(da$feature_ID), 5))
    cat("Tissues:", paste(unique(da$tissue), collapse = ", "), "\n")
    cat("Sexes:", paste(unique(da$sex), collapse = ", "), "\n")
    cat("Comparison groups:", paste(unique(da$comparison_group), collapse = ", "), "\n")
    cat("Total features:", n_distinct(da$feature_ID), "\n")
    cat("Significant (adj_p < 0.05):", sum(da$adj_p_value < 0.05, na.rm = TRUE), "rows\n")
  }
}

# Also check: can METHYL features be mapped to genes?
cat("\n=== METHYL feature → gene mapping check ===\n")
if (length(methyl_tables) > 0) {
  e <- new.env()
  data(list = methyl_tables[1], package = "MotrpacRatTraining6moData", envir = e)
  da <- get(methyl_tables[1], envir = e)

  methyl_features <- unique(da$feature_ID)
  cat("Total METHYL features in first table:", length(methyl_features), "\n")

  mapped <- FEATURE_TO_GENE %>%
    filter(feature_ID %in% methyl_features)
  cat("METHYL features found in FEATURE_TO_GENE:", n_distinct(mapped$feature_ID), "\n")
  cat("Unique genes mapped:", n_distinct(mapped$gene_symbol), "\n")

  if (nrow(mapped) > 0) {
    cat("Sample mappings:\n")
    print(head(mapped %>% select(feature_ID, gene_symbol), 10))
  }
}

# ----------------------------------------------------------
# 9. ATAC assay inventory (brief)
# ----------------------------------------------------------
cat("\n\n=== ATAC DA TABLE INVENTORY ===\n")
atac_tables <- da_tables[grepl("^ATAC_", da_tables)]
cat("ATAC DA tables:", length(atac_tables), "\n")
cat(paste(" ", atac_tables), sep = "\n")
cat("\n")

if (length(atac_tables) > 0) {
  e <- new.env()
  data(list = atac_tables[1], package = "MotrpacRatTraining6moData", envir = e)
  da <- get(atac_tables[1], envir = e)
  cat("First table dimensions:", nrow(da), "x", ncol(da), "\n")
  cat("Columns:", paste(names(da), collapse = ", "), "\n")
  cat("Sample feature_IDs:\n")
  print(head(unique(da$feature_ID), 5))

  atac_mapped <- FEATURE_TO_GENE %>%
    filter(feature_ID %in% unique(da$feature_ID))
  cat("ATAC features in FEATURE_TO_GENE:", n_distinct(atac_mapped$feature_ID), "\n")
  cat("Unique genes:", n_distinct(atac_mapped$gene_symbol), "\n")
}

# ----------------------------------------------------------
# 10. Save raw results
# ----------------------------------------------------------
cat("\n\nSaving data/aging_recon.csv...\n")
write_csv(results, "data/aging_recon.csv")
cat("Done:", nrow(results), "rows\n")
