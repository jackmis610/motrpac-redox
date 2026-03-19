# ============================================================
# 18: Txnrd2 vs Nnt cross-tissue analysis (TRNSCRPT, 8w, male/female)
# Key question: does Txnrd2 show inverse correlation with ETS
# like Nnt, suggesting NADPH-consuming enzymes scale inversely
# with ETS expansion?
# Produces: data/txnrd2_nnt_correlation.csv
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

source("scripts/neufer_geneset.R")

# ----------------------------------------------------------
# 1. Load neufer_da_all.csv
# ----------------------------------------------------------
cat("Loading data/neufer_da_all.csv...\n")
da_all <- read_csv("data/neufer_da_all.csv", show_col_types = FALSE)
cat("Total rows:", nrow(da_all), "\n\n")

# ----------------------------------------------------------
# 2. Define gene groups
# ----------------------------------------------------------
ets_genes <- c(unlist(neufer_geneset[c("ets_CI","ets_CIII","ets_CIV","cyt_c","atp_synthase")]))
focal_ets  <- c("Atp5f1b","Cox4i1","Ndufv1","Ndufs1","Uqcrc1")
buffering_genes <- c("Txnrd2","Nnt")

cat("ETS genes (", length(ets_genes), "): ", paste(focal_ets, collapse=", "), " + ...\n", sep="")
cat("Buffering genes:", paste(buffering_genes, collapse=", "), "\n\n")

# ----------------------------------------------------------
# Helper: compute tissue-level mean logFC for a gene × sex × time
# ----------------------------------------------------------
get_tissue_means <- function(data, genes, sex_sel, time_sel = "8w") {
  data %>%
    filter(
      assay == "TRNSCRPT",
      comparison_group == time_sel,
      sex == sex_sel,
      gene_symbol %in% genes
    ) %>%
    group_by(tissue, gene_symbol) %>%
    summarise(mean_logFC = mean(logFC, na.rm = TRUE), .groups = "drop")
}

# ----------------------------------------------------------
# 3. Build tissue-level means: males 8w
# ----------------------------------------------------------
cat("=== MALES, 8w ===\n")

ets_m <- get_tissue_means(da_all, ets_genes, "male")
buf_m <- get_tissue_means(da_all, buffering_genes, "male")
focal_m <- get_tissue_means(da_all, focal_ets, "male")

# ETS index per tissue: mean across all ETS genes
ets_index_m <- ets_m %>%
  group_by(tissue) %>%
  summarise(ets_index = mean(mean_logFC, na.rm = TRUE), .groups = "drop")

tissues_m <- ets_index_m$tissue
cat("Tissues with male ETS data:", paste(tissues_m, collapse = ", "), "\n")
cat("BAT present:", "BAT" %in% tissues_m, "\n\n")

# ----------------------------------------------------------
# 4. Pearson r: buffering genes vs ETS index (males)
# ----------------------------------------------------------
compute_cors <- function(buf_data, ets_idx, sex_label, exclude_bat = FALSE) {
  results <- list()
  tissues_use <- ets_idx$tissue
  if (exclude_bat) tissues_use <- tissues_use[tissues_use != "BAT"]

  ets_sub <- ets_idx %>% filter(tissue %in% tissues_use)

  for (g in unique(buf_data$gene_symbol)) {
    buf_sub <- buf_data %>%
      filter(gene_symbol == g, tissue %in% tissues_use) %>%
      left_join(ets_sub, by = "tissue") %>%
      filter(!is.na(mean_logFC), !is.na(ets_index))

    if (nrow(buf_sub) < 3) {
      results[[length(results)+1]] <- tibble(
        gene = g, sex = sex_label,
        exclude_BAT = exclude_bat,
        n_tissues = nrow(buf_sub),
        r = NA_real_, p_value = NA_real_
      )
      next
    }

    ct <- cor.test(buf_sub$mean_logFC, buf_sub$ets_index, method = "pearson")
    results[[length(results)+1]] <- tibble(
      gene = g, sex = sex_label,
      exclude_BAT = exclude_bat,
      n_tissues = nrow(buf_sub),
      r = ct$estimate,
      p_value = ct$p.value
    )
  }
  bind_rows(results)
}

# Males vs ETS index — with and without BAT
cor_m_with    <- compute_cors(buf_m, ets_index_m, "male", exclude_bat = FALSE)
cor_m_without <- compute_cors(buf_m, ets_index_m, "male", exclude_bat = TRUE)

cat("--- Males: Txnrd2 & Nnt vs full ETS index ---\n")
cat("\nWith BAT:\n")
print(cor_m_with)
cat("\nWithout BAT:\n")
print(cor_m_without)

# ----------------------------------------------------------
# 5. Focal ETS genes — Pearson r vs buffering genes (males)
# ----------------------------------------------------------
cat("\n--- Males: buffering genes vs individual focal ETS genes ---\n")

focal_cor_results <- list()

for (buf_g in buffering_genes) {
  buf_sub <- buf_m %>% filter(gene_symbol == buf_g)

  for (ets_g in focal_ets) {
    ets_sub <- focal_m %>% filter(gene_symbol == ets_g)

    # With BAT
    combined_with <- buf_sub %>%
      left_join(ets_sub, by = "tissue", suffix = c("_buf","_ets")) %>%
      filter(!is.na(mean_logFC_buf), !is.na(mean_logFC_ets))

    if (nrow(combined_with) >= 3) {
      ct <- cor.test(combined_with$mean_logFC_buf, combined_with$mean_logFC_ets, method = "pearson")
      r_with <- ct$estimate; p_with <- ct$p.value
    } else { r_with <- NA; p_with <- NA }

    # Without BAT
    combined_wout <- combined_with %>% filter(tissue != "BAT")
    if (nrow(combined_wout) >= 3) {
      ct2 <- cor.test(combined_wout$mean_logFC_buf, combined_wout$mean_logFC_ets, method = "pearson")
      r_wout <- ct2$estimate; p_wout <- ct2$p.value
    } else { r_wout <- NA; p_wout <- NA }

    focal_cor_results[[length(focal_cor_results)+1]] <- tibble(
      buffering_gene = buf_g,
      ets_gene       = ets_g,
      sex            = "male",
      n_with_BAT     = nrow(combined_with),
      r_with_BAT     = r_with,
      p_with_BAT     = p_with,
      n_without_BAT  = nrow(combined_wout),
      r_without_BAT  = r_wout,
      p_without_BAT  = p_wout
    )
  }
}

focal_cors <- bind_rows(focal_cor_results)
cat("\n")
print(focal_cors, n = 30)

# ----------------------------------------------------------
# 6. Females vs ETS index
# ----------------------------------------------------------
cat("\n=== FEMALES, 8w ===\n")

ets_f   <- get_tissue_means(da_all, ets_genes, "female")
buf_f   <- get_tissue_means(da_all, buffering_genes, "female")

ets_index_f <- ets_f %>%
  group_by(tissue) %>%
  summarise(ets_index = mean(mean_logFC, na.rm = TRUE), .groups = "drop")

tissues_f <- ets_index_f$tissue
cat("Tissues with female ETS data:", paste(tissues_f, collapse = ", "), "\n")

cor_f_with    <- compute_cors(buf_f, ets_index_f, "female", exclude_bat = FALSE)
cor_f_without <- compute_cors(buf_f, ets_index_f, "female", exclude_bat = TRUE)

cat("\nWith BAT:\n")
print(cor_f_with)
cat("\nWithout BAT:\n")
print(cor_f_without)

# ----------------------------------------------------------
# 7. Comparison table
# ----------------------------------------------------------
cat("\n=== COMPARISON TABLE: Txnrd2 vs Nnt (ETS index correlations) ===\n")

comparison_table <- bind_rows(
  cor_m_with    %>% mutate(bat_context = "with_BAT"),
  cor_m_without %>% mutate(bat_context = "without_BAT"),
  cor_f_with    %>% mutate(bat_context = "with_BAT"),
  cor_f_without %>% mutate(bat_context = "without_BAT")
) %>%
  select(gene, sex, bat_context, n_tissues, r, p_value) %>%
  arrange(gene, sex, bat_context)

print(comparison_table, n = 30)

# Check direction agreement
cat("\n--- Direction comparison (Txnrd2 vs Nnt) ---\n")
dir_check <- comparison_table %>%
  group_by(sex, bat_context) %>%
  summarise(
    txnrd2_r    = r[gene == "Txnrd2"],
    nnt_r       = r[gene == "Nnt"],
    same_direction = sign(r[gene == "Txnrd2"]) == sign(r[gene == "Nnt"]),
    .groups = "drop"
  )
print(dir_check)

# ----------------------------------------------------------
# 8. Tissue-level data for transparency
# ----------------------------------------------------------
cat("\n=== Tissue-level data: Txnrd2 vs ETS index (males 8w) ===\n")
txnrd2_tissue <- buf_m %>%
  filter(gene_symbol == "Txnrd2") %>%
  left_join(ets_index_m, by = "tissue") %>%
  arrange(ets_index)
print(txnrd2_tissue)

cat("\n=== Tissue-level data: Nnt vs ETS index (males 8w) ===\n")
nnt_tissue <- buf_m %>%
  filter(gene_symbol == "Nnt") %>%
  left_join(ets_index_m, by = "tissue") %>%
  arrange(ets_index)
print(nnt_tissue)

# ----------------------------------------------------------
# 9. Save results
# ----------------------------------------------------------
dir.create("data", showWarnings = FALSE)

# Main output: ETS index correlations (all sex × BAT combos)
write_csv(comparison_table, "data/txnrd2_nnt_correlation.csv")
cat("\nSaved: data/txnrd2_nnt_correlation.csv\n")

# Also save focal gene correlations
write_csv(focal_cors, "data/txnrd2_nnt_focal_ets_correlations.csv")
cat("Saved: data/txnrd2_nnt_focal_ets_correlations.csv\n")

# ----------------------------------------------------------
# 10. Clear summary with interpretation
# ----------------------------------------------------------
cat("\n")
cat("=============================================================\n")
cat("SUMMARY & INTERPRETATION\n")
cat("=============================================================\n")

cat("\nKey question: Does Txnrd2 show INVERSE correlation with ETS\n")
cat("expansion across tissues (like Nnt), suggesting NADPH-consuming\n")
cat("enzymes generally scale inversely with ETS?\n\n")

for (g in c("Txnrd2","Nnt")) {
  cat(sprintf("--- %s ---\n", g))
  for (sx in c("male","female")) {
    r_w  <- comparison_table %>% filter(gene==g, sex==sx, bat_context=="with_BAT") %>% pull(r)
    r_wo <- comparison_table %>% filter(gene==g, sex==sx, bat_context=="without_BAT") %>% pull(r)
    p_wo <- comparison_table %>% filter(gene==g, sex==sx, bat_context=="without_BAT") %>% pull(p_value)
    n_wo <- comparison_table %>% filter(gene==g, sex==sx, bat_context=="without_BAT") %>% pull(n_tissues)
    dir_lbl <- if (!is.na(r_wo) && r_wo < 0) "INVERSE (negative)" else if (!is.na(r_wo)) "POSITIVE" else "NA"
    cat(sprintf("  %s: r_with_BAT=%.3f  r_without_BAT=%.3f (p=%.4f, n=%d tissues) => %s\n",
                sx,
                ifelse(length(r_w)==0 || is.na(r_w), NA, r_w),
                ifelse(length(r_wo)==0 || is.na(r_wo), NA, r_wo),
                ifelse(length(p_wo)==0 || is.na(p_wo), NA, p_wo),
                ifelse(length(n_wo)==0 || is.na(n_wo), 0L, n_wo),
                dir_lbl))
  }
}

cat("\nNote: BAT is flagged as unreliable in the regression figure due to\n")
cat("multi-feature mapping artifacts (see MEMORY.md rigor notes).\n")
cat("r_without_BAT is the primary interpretation metric.\n")
cat("=============================================================\n")
