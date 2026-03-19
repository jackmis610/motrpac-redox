# ============================================================
# 17: Supplementary Table 1 — Neufer Framework Gene Set
# Produces: data/supplementary_table1_gene_set.csv
# ============================================================

suppressPackageStartupMessages({
  library(MotrpacRatTraining6moData)
  library(tidyverse)
})

source("scripts/neufer_geneset.R")

# ----------------------------------------------------------
# 1. Load annotation
# ----------------------------------------------------------
cat("Loading FEATURE_TO_GENE...\n")
data(FEATURE_TO_GENE)

# ----------------------------------------------------------
# 2. Define subcategory mapping
# ----------------------------------------------------------
subcategory_map <- c(
  ets_CI        = "ETS",
  ets_CIII      = "ETS",
  ets_CIV       = "ETS",
  cyt_c         = "ETS",
  atp_synthase  = "ETS",
  nnt           = "Redox buffering",   # also NADPH source — assign to Redox buffering per spec
  gsh           = "Redox buffering",
  trx           = "Redox buffering",
  sod           = "Redox buffering",
  nadph_other   = "NADPH source",
  q_pool        = "Q-pool feeder",
  pdh           = "PDH complex",
  beta_ox       = "Beta-oxidation",
  ucp           = "UCP",
  shuttles      = "Shuttle"
)

in_ets_index_cats      <- c("ets_CI", "ets_CIII", "ets_CIV", "cyt_c", "atp_synthase")
in_buffering_index_cats <- c("gsh", "trx", "nnt", "sod")

# ----------------------------------------------------------
# 3. Build gene × category data frame
#    Rules:
#    - Gpd2 appears in q_pool AND shuttles → keep in q_pool
#    - Dld appears in pdh AND beta_ox → keep in pdh
# ----------------------------------------------------------
gene_cat_raw <- imap_dfr(neufer_geneset, function(genes, cat) {
  tibble(gene_symbol = genes, functional_category = cat)
})

# Apply deduplication: first occurrence wins (q_pool listed before shuttles,
# pdh listed before beta_ox in the geneset)
gene_cat <- gene_cat_raw %>%
  distinct(gene_symbol, .keep_all = TRUE)

cat("Unique genes in gene set:", nrow(gene_cat), "\n")

# ----------------------------------------------------------
# 4. Get best Ensembl + Entrez IDs from FEATURE_TO_GENE
#    For TRNSCRPT features: ensembl_gene is the feature_ID itself
#    Use first mapping hit per gene (most common strategy)
# ----------------------------------------------------------
ftg_neufer <- FEATURE_TO_GENE %>%
  filter(!is.na(gene_symbol), gene_symbol %in% all_neufer_genes) %>%
  select(gene_symbol, ensembl_gene, entrez_gene) %>%
  distinct()

# Summarise: take first ensembl_gene and first entrez_gene per gene
gene_ids <- ftg_neufer %>%
  group_by(gene_symbol) %>%
  summarise(
    ensembl_id = first(na.omit(ensembl_gene)),
    entrez_id  = first(na.omit(entrez_gene)),
    .groups = "drop"
  )

cat("Genes with annotation:", nrow(gene_ids), "\n")
missing <- setdiff(all_neufer_genes, gene_ids$gene_symbol)
cat("Genes NOT in FEATURE_TO_GENE:", paste(missing, collapse = ", "), "\n\n")

# ----------------------------------------------------------
# 5. Build final table
# ----------------------------------------------------------
supp_table <- gene_cat %>%
  left_join(gene_ids, by = "gene_symbol") %>%
  mutate(
    subcategory       = subcategory_map[functional_category],
    human_ortholog    = toupper(gene_symbol),
    mapping_status    = ifelse(is.na(ensembl_id) & gene_symbol %in% missing,
                               "not_mapped", "mapped"),
    in_ets_index      = functional_category %in% in_ets_index_cats,
    in_buffering_index = functional_category %in% in_buffering_index_cats
  ) %>%
  select(
    gene_symbol,
    functional_category,
    subcategory,
    ensembl_id,
    entrez_id,
    human_ortholog,
    mapping_status,
    in_ets_index,
    in_buffering_index
  )

cat("=== Supplementary Table 1 — Gene Set ===\n")
cat("Total rows:", nrow(supp_table), "\n\n")

# ----------------------------------------------------------
# 6. Summary by category
# ----------------------------------------------------------
cat("Rows per functional_category:\n")
summary_tbl <- supp_table %>%
  group_by(functional_category, subcategory) %>%
  summarise(
    n_genes    = n(),
    n_mapped   = sum(mapping_status == "mapped"),
    n_not_mapped = sum(mapping_status == "not_mapped"),
    .groups = "drop"
  ) %>%
  arrange(subcategory, functional_category)

print(summary_tbl, n = 30)

cat("\nMapping status overall:\n")
print(table(supp_table$mapping_status))

cat("\nin_ets_index summary:\n")
print(table(supp_table$in_ets_index))

cat("\nin_buffering_index summary:\n")
print(table(supp_table$in_buffering_index))

# ----------------------------------------------------------
# 7. Save
# ----------------------------------------------------------
dir.create("data", showWarnings = FALSE)
write_csv(supp_table, "data/supplementary_table1_gene_set.csv")
cat("\nSaved: data/supplementary_table1_gene_set.csv\n")
cat("Rows:", nrow(supp_table), "\n")

# Show first few rows as sanity check
cat("\nFirst 10 rows:\n")
print(supp_table %>% head(10))
