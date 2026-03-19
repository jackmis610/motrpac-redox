# ============================================================
# 05: Mito dynamics analysis
#   (1) Map genes through FEATURE_TO_GENE, query all DA tables
#   (2) Tissue × category heatmap at 8w
#   (3) ETS vs. Ring 2 / Ring 3 category correlations
#   (4) Metabolomics extraction: lactate, pyruvate, NAD+, TCA
# ============================================================

suppressPackageStartupMessages({
  library(MotrpacRatTraining6moData)
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
})

dir.create("figures", showWarnings = FALSE)
dir.create("data",    showWarnings = FALSE)

source("scripts/mito_dynamics_geneset.R")   # loads mito_dynamics_geneset, all_mito_dynamics_genes

# ── shared aesthetics (mirrors 03_figures.R / 04_rigor.R) ────────────────────
SEX_COL   <- c(Male = "#2166AC", Female = "#D6604D")
SEX_SHAPE <- c(Male = 16L,       Female = 17L)
DIVPAL    <- "RdBu"

# Ring membership for ordering and annotation
RING1_CATS <- c("mct_transport","ldh","mito_pyruvate_carrier","pdh","pdk_pdp","lactate_signaling")
RING2_CATS <- c("fission","fusion","mitophagy","biogenesis")
RING3_CATS <- c("mito_upr","integrated_stress","ampk_axis","calcium","mito_derived_signals")
BRIDGE_CAT <- "ros_signaling"

CAT_ORDER  <- c(RING1_CATS, RING2_CATS, RING3_CATS, BRIDGE_CAT)

CAT_LABELS <- c(
  mct_transport          = "MCT transport",
  ldh                    = "LDH",
  mito_pyruvate_carrier  = "MPC",
  pdh                    = "PDH complex",
  pdk_pdp                = "PDK / PDP",
  lactate_signaling      = "Lactate signaling",
  fission                = "Fission",
  fusion                 = "Fusion",
  mitophagy              = "Mitophagy",
  biogenesis             = "Biogenesis",
  mito_upr               = "mtUPR",
  integrated_stress      = "ISR",
  ampk_axis              = "AMPK axis",
  calcium                = "Ca\u00b2\u207a / VDAC",
  mito_derived_signals   = "Mitokines",
  ros_signaling          = "ROS signaling (bridge)"
)

RING_FILL <- c(
  mct_transport="#FEE0D2", ldh="#FEE0D2", mito_pyruvate_carrier="#FEE0D2",
  pdh="#FEE0D2", pdk_pdp="#FEE0D2", lactate_signaling="#FEE0D2",
  fission="#DEEBF7", fusion="#DEEBF7", mitophagy="#DEEBF7", biogenesis="#DEEBF7",
  mito_upr="#E5F5E0", integrated_stress="#E5F5E0", ampk_axis="#E5F5E0",
  calcium="#E5F5E0", mito_derived_signals="#E5F5E0",
  ros_signaling="#F2F0FB"
)

TISSUE_ORDER <- c("SKM-GN","SKM-VL","HEART","BAT","WAT-SC",
                  "LIVER","KIDNEY","LUNG","ADRNL",
                  "BLOOD","CORTEX","HIPPOC","HYPOTH","VENACV",
                  "COLON","SMLINT","SPLEEN","OVARY","TESTES")

theme_pub <- theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text       = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank(),
    legend.key.size  = unit(0.4, "cm")
  )

# ══════════════════════════════════════════════════════════════════════════════
# PART 1 — Build annotation lookup and query all DA tables
# ══════════════════════════════════════════════════════════════════════════════
cat("Building annotation lookup...\n")
data(FEATURE_TO_GENE)

gene_category <- imap_dfr(mito_dynamics_geneset, function(genes, cat) {
  tibble(gene_symbol = genes, framework_category = cat)
}) %>% distinct(gene_symbol, .keep_all = TRUE)

mito_annotation <- FEATURE_TO_GENE %>%
  filter(!is.na(gene_symbol), gene_symbol %in% all_mito_dynamics_genes) %>%
  select(feature_ID, gene_symbol, ensembl_gene, entrez_gene) %>%
  distinct() %>%
  left_join(gene_category, by = "gene_symbol")

cat("Mito dynamics features in annotation:", nrow(mito_annotation), "\n")
cat("Unique genes mapped:", n_distinct(mito_annotation$gene_symbol),
    "of", length(all_mito_dynamics_genes), "\n")

missing <- setdiff(all_mito_dynamics_genes, mito_annotation$gene_symbol)
if (length(missing)) cat("Genes not found:", paste(missing, collapse=", "), "\n\n")

# Discover DA tables (skip METAB, IMMUNO, ATAC, METHYL)
all_items  <- data(package = "MotrpacRatTraining6moData")$results[, "Item"]
da_tables  <- all_items[grepl("_DA$", all_items)]
da_tables  <- da_tables[!grepl("^METAB_|^IMMUNO_|^ATAC_|^METHYL_", da_tables)]

cat("Querying", length(da_tables), "DA tables...\n")

query_da <- function(table_name) {
  e <- new.env()
  data(list = table_name, package = "MotrpacRatTraining6moData", envir = e)
  da <- get(table_name, envir = e)
  hits <- da %>%
    inner_join(mito_annotation, by = "feature_ID") %>%
    select(
      gene_symbol, framework_category, feature_ID, ensembl_gene,
      assay, tissue, sex, comparison_group,
      logFC, adj_p_value, p_value,
      any_of(c("shrunk_logFC", "zscore", "selection_fdr", "tscore"))
    )
  if (nrow(hits) > 0)
    cat(sprintf("  %-30s %d rows\n", table_name, nrow(hits)))
  hits
}

results_list <- map(da_tables, safely(query_da))

mito_da <- map(results_list, "result") %>% compact() %>% bind_rows()

walk2(da_tables, results_list, function(nm, res) {
  if (!is.null(res$error))
    cat("  ERROR in", nm, ":", conditionMessage(res$error), "\n")
})

cat(sprintf("\nTotal: %d rows, %d genes\n",
            nrow(mito_da), n_distinct(mito_da$gene_symbol)))

sig_da <- mito_da %>% filter(adj_p_value < 0.05)
cat("Significant (adj_p < 0.05):", nrow(sig_da), "rows,",
    n_distinct(sig_da$gene_symbol), "genes\n\n")

write_csv(mito_da,  "data/mito_dynamics_da_all.csv")
write_csv(sig_da,   "data/mito_dynamics_da_significant.csv")
cat("Saved data/mito_dynamics_da_all.csv\n")
cat("Saved data/mito_dynamics_da_significant.csv\n\n")

# Significant summary
sig_summary <- sig_da %>%
  group_by(framework_category, gene_symbol, assay, tissue, sex) %>%
  summarise(n_sig_tp = n(), min_adj_p = min(adj_p_value),
            max_abs_logFC = max(abs(logFC)), .groups = "drop") %>%
  arrange(framework_category, min_adj_p)
cat("=== Top significant hits ===\n")
print(sig_summary, n = 30)

# ══════════════════════════════════════════════════════════════════════════════
# PART 2 — Tissue × category heatmap (TRNSCRPT, 8w)
# ══════════════════════════════════════════════════════════════════════════════
cat("\nBuilding heatmap...\n")

heat_dat <- mito_da %>%
  filter(assay == "TRNSCRPT", comparison_group == "8w") %>%
  mutate(
    framework_category = factor(framework_category, levels = CAT_ORDER),
    tissue             = factor(tissue,             levels = TISSUE_ORDER),
    sex_label          = str_to_title(sex)
  ) %>%
  group_by(tissue, framework_category, sex_label) %>%
  summarise(
    mean_logFC = mean(logFC, na.rm = TRUE),
    n_genes    = n_distinct(gene_symbol),
    n_sig      = sum(adj_p_value < 0.05, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(sig_label = case_when(
    n_sig >= 3 ~ "***", n_sig == 2 ~ "**", n_sig == 1 ~ "*", TRUE ~ ""
  ))

lim <- ceiling(max(abs(heat_dat$mean_logFC), na.rm = TRUE) * 10) / 10

fig_heat <- ggplot(heat_dat,
                   aes(x = framework_category, y = tissue, fill = mean_logFC)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sig_label), size = 2.6, vjust = 0.8, colour = "grey20") +
  scale_fill_distiller(palette = DIVPAL, limits = c(-lim, lim),
                       direction = -1, na.value = "grey88",
                       name = "Mean log\u2082FC\nvs. sedentary") +
  scale_x_discrete(labels = CAT_LABELS) +
  scale_y_discrete(limits = rev(levels(heat_dat$tissue))) +
  # Ring separator lines
  geom_vline(xintercept = c(6.5, 10.5, 15.5),
             colour = "grey30", linewidth = 0.6, linetype = "solid") +
  annotate("text", x = 3.5, y = 20.2, label = "Ring 1: Lactate axis",
           size = 2.8, fontface = "bold", hjust = 0.5, vjust = 0) +
  annotate("text", x = 8.5, y = 20.2, label = "Ring 2: Dynamics",
           size = 2.8, fontface = "bold", hjust = 0.5, vjust = 0) +
  annotate("text", x = 13,  y = 20.2, label = "Ring 3: Signaling",
           size = 2.8, fontface = "bold", hjust = 0.5, vjust = 0) +
  annotate("text", x = 16,  y = 20.2, label = "Bridge",
           size = 2.8, fontface = "bold", hjust = 0.5, vjust = 0) +
  facet_wrap(~ sex_label, ncol = 1) +
  coord_cartesian(clip = "off") +
  labs(
    title    = "Mito Dynamics Gene Set — 8-Week Training Response",
    subtitle = "Mean log\u2082FC (TRNSCRPT) vs. sedentary; * \u22651 gene adj_p < 0.05",
    x = NULL, y = NULL
  ) +
  theme_pub +
  theme(
    axis.text.x     = element_text(angle = 40, hjust = 1, size = 8),
    axis.text.y     = element_text(size = 8.5),
    panel.border    = element_blank(),
    legend.position = "right",
    plot.margin     = margin(t = 22, r = 10, b = 5, l = 5)
  )

ggsave("figures/fig_mito_heatmap.pdf",
       fig_heat, width = 13, height = 11, device = cairo_pdf)
ggsave("figures/fig_mito_heatmap.png",
       fig_heat, width = 13, height = 11, dpi = 200)
cat("Saved figures/fig_mito_heatmap\n")

# ══════════════════════════════════════════════════════════════════════════════
# PART 3 — ETS vs. Ring 2 / Ring 3 category correlations across tissues
# ══════════════════════════════════════════════════════════════════════════════
cat("\nRunning ETS vs. Ring 2/3 correlations...\n")

# Load neufer ETS index (already computed; recompute inline for independence)
neufer_da <- read_csv("data/neufer_da_all.csv", show_col_types = FALSE)

ets_index <- neufer_da %>%
  filter(assay == "TRNSCRPT", comparison_group == "8w",
         framework_category %in% c("ets_CI","ets_CIII","ets_CIV")) %>%
  group_by(tissue, sex) %>%
  summarise(ets_mean = mean(logFC, na.rm = TRUE), .groups = "drop")

# Ring 2 and Ring 3 category means per tissue × sex (TRNSCRPT, 8w)
ring23_means <- mito_da %>%
  filter(assay == "TRNSCRPT", comparison_group == "8w",
         framework_category %in% c(RING2_CATS, RING3_CATS)) %>%
  group_by(framework_category, tissue, sex) %>%
  summarise(cat_mean = mean(logFC, na.rm = TRUE), .groups = "drop")

ring23_ets <- ring23_means %>%
  left_join(ets_index, by = c("tissue","sex")) %>%
  filter(!is.na(ets_mean)) %>%
  mutate(
    sex_label          = str_to_title(sex),
    framework_category = factor(framework_category,
                                levels = c(RING2_CATS, RING3_CATS))
  )

# Pearson r per category × sex
ring23_cor <- ring23_ets %>%
  group_by(framework_category, sex_label) %>%
  summarise(
    r    = cor(ets_mean, cat_mean, use = "complete.obs"),
    pval = tryCatch(cor.test(ets_mean, cat_mean)$p.value, error = function(e) NA_real_),
    n    = sum(!is.na(cat_mean) & !is.na(ets_mean)),
    .groups = "drop"
  ) %>%
  mutate(
    sig   = case_when(pval < 0.001 ~ "***", pval < 0.01 ~ "**",
                      pval < 0.05 ~ "*", TRUE ~ ""),
    label = sprintf("r=%.2f%s", r, sig)
  )

cat("\nETS vs Ring 2/3 category correlations (TRNSCRPT, 8w, n=19 tissues):\n")
print(ring23_cor %>% arrange(sex_label, desc(abs(r))), n = 40)

write_csv(ring23_cor, "data/mito_ets_ring23_correlations.csv")
cat("Saved data/mito_ets_ring23_correlations.csv\n")

# Plot: lollipop of r values per category, faceted by sex
fig_cor <- ggplot(ring23_cor,
                  aes(x = r, y = reorder(framework_category, r),
                      colour = sex_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60",
             linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = r,
                   y = framework_category, yend = framework_category),
               linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 3.5) +
  geom_text(aes(label = label),
            nudge_x = ifelse(ring23_cor$r > 0, 0.04, -0.04),
            hjust    = ifelse(ring23_cor$r > 0, 0, 1),
            size = 2.8, show.legend = FALSE) +
  scale_colour_manual(values = SEX_COL, name = "Sex") +
  scale_y_discrete(labels = CAT_LABELS) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
  facet_wrap(~ sex_label, ncol = 2) +
  labs(
    title    = "ETS Expansion vs. Ring 2 / Ring 3 Category Expression",
    subtitle = "Pearson r across 19 tissues (TRNSCRPT, 8w). * p<0.05 ** p<0.01 *** p<0.001",
    x        = "Pearson r  (ETS mean log\u2082FC ~ category mean log\u2082FC)",
    y        = NULL
  ) +
  theme_pub +
  theme(legend.position = "none",
        panel.grid.major.y = element_line(colour = "grey93"))

ggsave("figures/fig_mito_ets_ring23_cor.pdf",
       fig_cor, width = 11, height = 6, device = cairo_pdf)
ggsave("figures/fig_mito_ets_ring23_cor.png",
       fig_cor, width = 11, height = 6, dpi = 200)
cat("Saved figures/fig_mito_ets_ring23_cor\n")

# Scatter for top correlated categories
top_cats <- ring23_cor %>%
  filter(pval < 0.05) %>%
  pull(framework_category) %>%
  unique()

if (length(top_cats) > 0) {
  scatter_dat <- ring23_ets %>%
    filter(framework_category %in% top_cats)

  cor_ann <- ring23_cor %>%
    filter(framework_category %in% top_cats) %>%
    mutate(framework_category = factor(framework_category,
                                       levels = c(RING2_CATS, RING3_CATS)))

  fig_cor_scatter <- ggplot(scatter_dat,
                            aes(x = ets_mean, y = cat_mean,
                                colour = sex_label, shape = sex_label,
                                label = tissue)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70",
               linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70",
               linewidth = 0.3) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.15) +
    geom_point(size = 2.5) +
    geom_text_repel(size = 2.4, max.overlaps = 15, seed = 42,
                    segment.size = 0.25, show.legend = FALSE) +
    geom_text(
      data = cor_ann,
      aes(x = -Inf, y = Inf, label = label, colour = sex_label),
      hjust = -0.1, vjust = 1.5, size = 2.8, inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    scale_colour_manual(values = SEX_COL,   name = "Sex") +
    scale_shape_manual( values = SEX_SHAPE, name = "Sex") +
    facet_grid(sex_label ~ framework_category,
               labeller = labeller(framework_category = CAT_LABELS)) +
    labs(
      title    = "ETS Expansion vs. Significantly Correlated Ring 2/3 Categories",
      subtitle = "Each point = one tissue; sig. categories only (p < 0.05 in \u22651 sex)",
      x        = "ETS mean log\u2082FC (CI+III+IV, TRNSCRPT 8w)",
      y        = "Category mean log\u2082FC"
    ) +
    theme_pub +
    theme(legend.position = "bottom",
          strip.text.x    = element_text(size = 7.5))

  ggsave("figures/fig_mito_ets_scatter.pdf",
         fig_cor_scatter, width = 11, height = 7, device = cairo_pdf)
  ggsave("figures/fig_mito_ets_scatter.png",
         fig_cor_scatter, width = 11, height = 7, dpi = 200)
  cat("Saved figures/fig_mito_ets_scatter\n")
} else {
  cat("No Ring 2/3 categories significantly correlated with ETS at p < 0.05\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# PART 4 — Metabolomics extraction: lactate, pyruvate, NAD+, TCA across tissues
# ══════════════════════════════════════════════════════════════════════════════
cat("\nExtracting metabolomics data...\n")

# Canonical metabolite groups (case-insensitive matching, then harmonize name)
metab_groups <- tribble(
  ~feature_ID,             ~canonical,          ~metab_group,
  "lactate",               "Lactate",            "Lactate/Pyruvate",
  "Lactate",               "Lactate",            "Lactate/Pyruvate",
  "pyruvate",              "Pyruvate",           "Lactate/Pyruvate",
  "Pyruvate",              "Pyruvate",           "Lactate/Pyruvate",
  "NAD+",                  "NAD+",               "NAD redox",
  "NADH",                  "NADH",               "NAD redox",
  "NADP+",                 "NADP+",              "NAD redox",
  "NADPH",                 "NADPH",              "NAD redox",
  "FAD",                   "FAD",                "NAD redox",
  "succinate",             "Succinate",          "TCA cycle",
  "Succinate",             "Succinate",          "TCA cycle",
  "fumarate",              "Fumarate",           "TCA cycle",
  "Fumarate",              "Fumarate",           "TCA cycle",
  "malate",                "Malate",             "TCA cycle",
  "Malate",                "Malate",             "TCA cycle",
  "citrate",               "Citrate",            "TCA cycle",
  "Citrate",               "Citrate",            "TCA cycle",
  "citrate+isocitrate",    "Citrate+Isocitrate", "TCA cycle",
  "Isocitrate",            "Isocitrate",         "TCA cycle",
  "isocitrate",            "Isocitrate",         "TCA cycle",
  "alpha-ketoglutarate",   "alpha-KG",           "TCA cycle",
  "2-oxoglutarate",        "alpha-KG",           "TCA cycle",
  "phosphoenolpyruvate",   "PEP",                "Glycolysis",
  "glucose",               "Glucose",            "Glycolysis",
  "Glucose",               "Glucose",            "Glycolysis",
  "glucose-6-phosphate",   "G6P",                "Glycolysis",
  "fructose-6-phosphate",  "F6P",                "Glycolysis"
)

metab_group_order <- c("Lactate/Pyruvate","NAD redox","TCA cycle","Glycolysis")

# Load all METAB DA tables and filter
metab_tables <- all_items[grepl("^METAB_.*_DA$", all_items)]

cat("Extracting from", length(metab_tables), "METAB tables...\n")
metab_list <- map(metab_tables, function(tbl) {
  e <- new.env()
  data(list = tbl, package = "MotrpacRatTraining6moData", envir = e)
  da <- get(tbl, envir = e)
  hits <- da %>%
    inner_join(metab_groups, by = "feature_ID") %>%
    select(canonical, metab_group, feature_ID, tissue, sex,
           comparison_group, logFC, adj_p_value, p_value,
           any_of(c("reference_average_intensity","comparison_average_intensity")))
  if (nrow(hits) > 0)
    cat(sprintf("  %-22s %d rows (%s)\n", tbl, nrow(hits), unique(da$tissue)))
  hits
})

metab_all <- bind_rows(metab_list) %>%
  # Average duplicate platform entries for same metabolite × tissue × sex × timepoint
  group_by(canonical, metab_group, tissue, sex, comparison_group) %>%
  summarise(
    logFC           = mean(logFC,       na.rm = TRUE),
    adj_p_value     = min(adj_p_value,  na.rm = TRUE),
    p_value         = min(p_value,      na.rm = TRUE),
    n_platform_entries = n(),
    .groups = "drop"
  ) %>%
  mutate(
    sex_label          = str_to_title(sex),
    comparison_group   = factor(comparison_group, levels = c("1w","2w","4w","8w")),
    metab_group        = factor(metab_group,       levels = metab_group_order),
    tissue             = factor(tissue,            levels = TISSUE_ORDER),
    sig_dot            = if_else(adj_p_value < 0.05, "\u25cf", "")
  )

cat("\nTotal metabolomics rows extracted:", nrow(metab_all), "\n")
cat("Tissues:", paste(sort(unique(metab_all$tissue)), collapse=", "), "\n")
cat("Metabolites:", paste(sort(unique(metab_all$canonical)), collapse=", "), "\n")

write_csv(metab_all, "data/tca_nad_lactate_metab_all.csv")
cat("Saved data/tca_nad_lactate_metab_all.csv\n\n")

# Significant metabolite hits
metab_sig <- metab_all %>% filter(adj_p_value < 0.05)
cat("Significant metabolite hits (adj_p < 0.05):\n")
print(metab_sig %>% arrange(adj_p_value) %>%
      select(canonical, tissue, sex_label, comparison_group, logFC, adj_p_value),
      n = 40)
write_csv(metab_sig, "data/tca_nad_lactate_metab_significant.csv")
cat("Saved data/tca_nad_lactate_metab_significant.csv\n")

# ── Metabolomics visualization ────────────────────────────────────────────────
# Heatmap: tissue × metabolite at 8w (mean across sex, or split)
metab_8w <- metab_all %>%
  filter(comparison_group == "8w") %>%
  mutate(sig_dot = if_else(adj_p_value < 0.05, "\u25cf", ""))

metab_lim <- ceiling(max(abs(metab_8w$logFC), na.rm = TRUE) * 10) / 10

fig_metab_heat <- ggplot(metab_8w,
                         aes(x = canonical, y = tissue, fill = logFC)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = sig_dot), size = 2.4, colour = "grey20", vjust = 0.75) +
  scale_fill_distiller(
    palette = DIVPAL, limits = c(-metab_lim, metab_lim),
    direction = -1, na.value = "grey88",
    name = "log\u2082FC\nvs. sedentary"
  ) +
  scale_y_discrete(limits = rev(intersect(levels(metab_8w$tissue),
                                          unique(metab_8w$tissue)))) +
  facet_grid(sex_label ~ metab_group, scales = "free_x", space = "free_x") +
  labs(
    title    = "Lactate / Pyruvate / NAD+ / TCA Metabolites — 8-Week Response",
    subtitle = "\u25cf adj_p < 0.05 (MoTrPAC FDR)",
    x = NULL, y = NULL
  ) +
  theme_pub +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7.5),
    axis.text.y     = element_text(size = 8),
    panel.border    = element_rect(colour = "grey80", fill = NA),
    legend.position = "right"
  )

ggsave("figures/fig_metab_heatmap_8w.pdf",
       fig_metab_heat, width = 14, height = 8, device = cairo_pdf)
ggsave("figures/fig_metab_heatmap_8w.png",
       fig_metab_heat, width = 14, height = 8, dpi = 200)
cat("\nSaved figures/fig_metab_heatmap_8w\n")

# Line trajectories for tissues with full coverage
full_tissues <- metab_all %>%
  group_by(tissue, sex_label) %>%
  summarise(n_metabs = n_distinct(canonical), .groups = "drop") %>%
  filter(n_metabs >= 10) %>%
  pull(tissue) %>%
  unique()

cat("Tissues with full metabolomics coverage (>=10 metabolites):",
    paste(full_tissues, collapse=", "), "\n")

metab_traj <- metab_all %>%
  filter(tissue %in% full_tissues,
         canonical %in% c("Lactate","Pyruvate","NAD+","NADH","NADP+","NADPH",
                          "Succinate","Malate","Citrate","Fumarate")) %>%
  mutate(canonical = factor(canonical,
                            levels = c("Lactate","Pyruvate",
                                       "NAD+","NADH","NADP+","NADPH",
                                       "Succinate","Malate","Citrate","Fumarate")))

fig_metab_traj <- ggplot(metab_traj,
                         aes(x = comparison_group, y = logFC,
                             group = interaction(tissue, sex_label),
                             colour = tissue, linetype = sex_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70",
             linewidth = 0.3) +
  geom_line(linewidth = 0.65, alpha = 0.8) +
  geom_point(aes(shape = sex_label), size = 1.8) +
  geom_text(aes(label = sig_dot), nudge_x = 0.18, size = 2.8,
            show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2", name = "Tissue") +
  scale_linetype_manual(values = c(Male = "solid", Female = "dashed"),
                        name = "Sex") +
  scale_shape_manual(values  = c(Male = 16L, Female = 17L), name = "Sex") +
  facet_wrap(~ canonical, scales = "free_y", ncol = 5) +
  labs(
    title    = "Key Metabolite Trajectories Across Training (tissues with full coverage)",
    subtitle = "\u25cf adj_p < 0.05; solid = male, dashed = female",
    x = "Training duration", y = "log\u2082FC vs. sedentary"
  ) +
  theme_pub +
  theme(legend.position = "bottom",
        legend.box       = "horizontal")

ggsave("figures/fig_metab_trajectories.pdf",
       fig_metab_traj, width = 14, height = 7, device = cairo_pdf)
ggsave("figures/fig_metab_trajectories.png",
       fig_metab_traj, width = 14, height = 7, dpi = 200)
cat("Saved figures/fig_metab_trajectories\n")

cat("\nDone. Summary:\n")
cat("  data/mito_dynamics_da_all.csv\n")
cat("  data/mito_dynamics_da_significant.csv\n")
cat("  data/mito_ets_ring23_correlations.csv\n")
cat("  data/tca_nad_lactate_metab_all.csv\n")
cat("  data/tca_nad_lactate_metab_significant.csv\n")
cat("  figures/fig_mito_heatmap\n")
cat("  figures/fig_mito_ets_ring23_cor\n")
cat("  figures/fig_mito_ets_scatter  (if sig categories found)\n")
cat("  figures/fig_metab_heatmap_8w\n")
cat("  figures/fig_metab_trajectories\n")
