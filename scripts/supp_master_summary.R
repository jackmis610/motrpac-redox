# ============================================================
# 07: Master tissue profile table — Supplementary Table 1
#
# One row per tissue × sex (then averaged for the single-row view):
#   ETS mean logFC, buffering mean logFC, buffering/ETS ratio,
#   top Ring2/3 category, top significant gene, trx-prot concordance
#
# Outputs:
#   data/tissue_profiles_master.csv   — full wide table (sex-split)
#   data/tissue_profiles_master_avg.csv — sex-averaged summary
#   figures/fig_supp_table1_*         — formatted visual table
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(gridExtra)
})

dir.create("data",    showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ── Constants ─────────────────────────────────────────────────────────────────
TISSUE_ORDER <- c("SKM-GN","SKM-VL","HEART","BAT","WAT-SC",
                  "LIVER","KIDNEY","LUNG","ADRNL","BLOOD",
                  "CORTEX","HIPPOC","HYPOTH","VENACV",
                  "COLON","SMLINT","SPLEEN","OVARY","TESTES")

# Tissue code → display name lookup
TC2NAME <- c(SKMGN="SKM-GN", SKMVL="SKM-VL", HEART="HEART",  BAT="BAT",
             WATSC="WAT-SC", LIVER="LIVER",  KIDNEY="KIDNEY", LUNG="LUNG",
             ADRNL="ADRNL",  BLOOD="BLOOD",  CORTEX="CORTEX", HIPPOC="HIPPOC",
             HYPOTH="HYPOTH", VENACV="VENACV", COLON="COLON",
             SMLINT="SMLINT", SPLEEN="SPLEEN", OVARY="OVARY", TESTES="TESTES")

# Assay coverage per tissue (from earlier audit)
PROT_TISSUES    <- c("SKM-GN","HEART","KIDNEY","LIVER","LUNG","CORTEX","WAT-SC")
PHOSPHO_TISSUES <- c("SKM-GN","HEART","KIDNEY","LIVER","LUNG","CORTEX","WAT-SC")
ACETYL_TISSUES  <- c("HEART","LIVER")
UBIQ_TISSUES    <- c("HEART","LIVER")
METAB_TISSUES   <- c("SKM-GN","BAT","WAT-SC","LIVER","KIDNEY","LUNG","HEART",
                      "HIPPOC","ADRNL","COLON","CORTEX","HYPOTH","SMLINT",
                      "SPLEEN","OVARY","TESTES","VENACV","SKM-VL")

ETS_CATS <- c("ets_CI","ets_CIII","ets_CIV")
BUF_CATS <- c("gsh","trx","nnt","sod")
R23_CATS <- c("fission","fusion","mitophagy","biogenesis",
              "mito_upr","integrated_stress","ampk_axis","calcium",
              "mito_derived_signals")
CAT_LABELS <- c(
  fission="Fission", fusion="Fusion", mitophagy="Mitophagy",
  biogenesis="Biogenesis", mito_upr="mtUPR", integrated_stress="ISR",
  ampk_axis="AMPK", calcium="Ca²⁺/VDAC", mito_derived_signals="Mitokines",
  ros_signaling="ROS bridge"
)

theme_pub <- theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill="grey92", colour=NA),
        strip.text = element_text(face="bold", size=8))

# ── Load data ─────────────────────────────────────────────────────────────────
cat("Loading data...\n")
neufer  <- read_csv("data/neufer_da_all.csv",        show_col_types = FALSE)
mito    <- read_csv("data/mito_dynamics_da_all.csv", show_col_types = FALSE)
concord <- read_csv("data/trx_prot_concordance.csv", show_col_types = FALSE)

all_genes <- bind_rows(neufer, mito) %>%
  mutate(sex_label = str_to_title(sex))

# ══════════════════════════════════════════════════════════════════════════════
# BLOCK 1 — ETS and buffering indices (TRNSCRPT, 8w)
# ══════════════════════════════════════════════════════════════════════════════
cat("Computing ETS/buffering indices...\n")

trx8 <- all_genes %>%
  filter(assay == "TRNSCRPT", comparison_group == "8w")

ets_stats <- trx8 %>%
  filter(framework_category %in% ETS_CATS) %>%
  group_by(tissue, sex_label) %>%
  summarise(
    ets_mean    = mean(logFC, na.rm = TRUE),
    ets_sd      = sd(logFC,   na.rm = TRUE),
    ets_n_genes = n_distinct(gene_symbol),
    ets_n_sig   = sum(adj_p_value < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

buf_stats <- trx8 %>%
  filter(framework_category %in% BUF_CATS) %>%
  group_by(tissue, sex_label) %>%
  summarise(
    buf_mean    = mean(logFC, na.rm = TRUE),
    buf_sd      = sd(logFC,   na.rm = TRUE),
    buf_n_genes = n_distinct(gene_symbol),
    buf_n_sig   = sum(adj_p_value < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

# ══════════════════════════════════════════════════════════════════════════════
# BLOCK 2 — Ring 2/3 top category per tissue (highest absolute mean logFC)
# ══════════════════════════════════════════════════════════════════════════════
cat("Finding top Ring 2/3 category per tissue...\n")

r23_per_tissue <- trx8 %>%
  filter(framework_category %in% R23_CATS) %>%
  group_by(tissue, sex_label, framework_category) %>%
  summarise(cat_mean = mean(logFC, na.rm = TRUE),
            cat_n_sig = sum(adj_p_value < 0.05, na.rm = TRUE),
            .groups = "drop")

# Top category = highest absolute mean logFC (averaged across sex first)
top_r23 <- r23_per_tissue %>%
  group_by(tissue, framework_category) %>%
  summarise(cat_mean_avg = mean(cat_mean, na.rm = TRUE),
            cat_n_sig_total = sum(cat_n_sig), .groups = "drop") %>%
  group_by(tissue) %>%
  slice_max(abs(cat_mean_avg), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(tissue,
            top_r23_cat   = CAT_LABELS[framework_category],
            top_r23_mean  = round(cat_mean_avg, 3),
            top_r23_n_sig = cat_n_sig_total)

# ══════════════════════════════════════════════════════════════════════════════
# BLOCK 3 — Top significant gene per tissue (lowest adj_p, 8w TRNSCRPT,
#           across all Neufer + mito categories; prefer genes with adj_p < 0.05)
# ══════════════════════════════════════════════════════════════════════════════
cat("Finding top significant gene per tissue...\n")

top_gene <- trx8 %>%
  group_by(tissue) %>%
  arrange(adj_p_value, p_value) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(tissue,
            top_gene         = gene_symbol,
            top_gene_cat     = framework_category,
            top_gene_lfc     = round(logFC, 3),
            top_gene_adj_p   = signif(adj_p_value, 3),
            top_gene_sex     = sex_label,
            top_gene_dataset = if_else(gene_symbol %in% unique(neufer$gene_symbol),
                                       "Neufer", "MitoDyn"))

# ══════════════════════════════════════════════════════════════════════════════
# BLOCK 4 — Transcript-protein concordance (ETS, 8w; 7 tissues)
# ══════════════════════════════════════════════════════════════════════════════
cat("Computing concordance averages...\n")

concord_avg <- concord %>%
  group_by(tissue) %>%
  summarise(
    trx_prot_r_all = round(mean(r_all, na.rm = TRUE), 2),
    trx_prot_r_ets = round(mean(r_ets, na.rm = TRUE), 2),
    trx_prot_ets_sig = sum(p_ets < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

# ══════════════════════════════════════════════════════════════════════════════
# BLOCK 5 — Assemble sex-split table
# ══════════════════════════════════════════════════════════════════════════════
cat("Assembling master table...\n")

sex_table <- ets_stats %>%
  left_join(buf_stats, by = c("tissue","sex_label")) %>%
  mutate(
    buf_ets_ratio = case_when(
      abs(ets_mean) < 0.01 ~ NA_real_,   # undefined when ETS ≈ 0
      TRUE                 ~ round(buf_mean / ets_mean, 3)
    ),
    ets_mean  = round(ets_mean, 3),
    buf_mean  = round(buf_mean, 3)
  )

# ══════════════════════════════════════════════════════════════════════════════
# BLOCK 6 — Sex-averaged master table (one row per tissue)
# ══════════════════════════════════════════════════════════════════════════════
master_wide <- sex_table %>%
  select(tissue, sex_label, ets_mean, buf_mean, buf_ets_ratio,
         ets_n_sig, buf_n_sig, ets_n_genes, buf_n_genes) %>%
  pivot_wider(
    names_from  = sex_label,
    values_from = c(ets_mean, buf_mean, buf_ets_ratio, ets_n_sig, buf_n_sig)
  ) %>%
  # sex-average for primary display columns
  mutate(
    ets_mean_avg = round((ets_mean_Male + ets_mean_Female) / 2, 3),
    buf_mean_avg = round((buf_mean_Male + buf_mean_Female) / 2, 3),
    buf_ets_ratio_male   = buf_ets_ratio_Male,
    buf_ets_ratio_female = buf_ets_ratio_Female,
    ets_n_sig_male   = ets_n_sig_Male,
    ets_n_sig_female = ets_n_sig_Female,
    buf_n_sig_male   = buf_n_sig_Male,
    buf_n_sig_female = buf_n_sig_Female,
    sex_diff_flag = abs(ets_mean_Male - ets_mean_Female) > 0.15
  ) %>%
  left_join(top_r23,      by = "tissue") %>%
  left_join(top_gene,     by = "tissue") %>%
  left_join(concord_avg,  by = "tissue") %>%
  # Data availability flags
  mutate(
    has_PROT    = tissue %in% PROT_TISSUES,
    has_PHOSPHO = tissue %in% PHOSPHO_TISSUES,
    has_ACETYL  = tissue %in% ACETYL_TISSUES,
    has_UBIQ    = tissue %in% UBIQ_TISSUES,
    has_METAB   = tissue %in% METAB_TISSUES,
    assay_coverage = paste0(
      if_else(has_PROT,    "P",  ""),
      if_else(has_PHOSPHO, "Ph", ""),
      if_else(has_ACETYL,  "A",  ""),
      if_else(has_UBIQ,    "U",  ""),
      if_else(has_METAB,   "M",  "")
    ),
    tissue = factor(tissue, levels = TISSUE_ORDER)
  ) %>%
  arrange(tissue)

# Clean column selection for CSV output
master_csv <- master_wide %>%
  select(
    tissue,
    ets_mean_male   = ets_mean_Male,
    ets_mean_female = ets_mean_Female,
    ets_mean_avg,
    ets_n_sig_male, ets_n_sig_female,
    buf_mean_male   = buf_mean_Male,
    buf_mean_female = buf_mean_Female,
    buf_mean_avg,
    buf_n_sig_male, buf_n_sig_female,
    buf_ets_ratio_male, buf_ets_ratio_female,
    top_r23_cat, top_r23_mean, top_r23_n_sig,
    top_gene, top_gene_cat, top_gene_lfc, top_gene_adj_p,
    top_gene_sex, top_gene_dataset,
    trx_prot_r_all, trx_prot_r_ets, trx_prot_ets_sig,
    assay_coverage,
    has_PROT, has_PHOSPHO, has_ACETYL, has_UBIQ, has_METAB,
    sex_diff_flag
  )

write_csv(master_csv, "data/tissue_profiles_master.csv")
cat("Saved data/tissue_profiles_master.csv\n")
cat("Rows:", nrow(master_csv), " | Cols:", ncol(master_csv), "\n\n")

# Preview
cat("=== Master table preview ===\n")
master_csv %>%
  select(tissue, ets_mean_avg, buf_mean_avg, buf_ets_ratio_male,
         top_r23_cat, top_gene, top_gene_adj_p,
         trx_prot_r_ets, assay_coverage) %>%
  as.data.frame() %>%
  print()

# ══════════════════════════════════════════════════════════════════════════════
# BLOCK 7 — Formatted visual table (Supplementary Figure)
# Three panels:
#   A) Heatmap columns: ETS mean (M/F), buffering mean (M/F), ratio bar
#   B) Annotation columns: top gene, top R23 category, concordance
#   C) Coverage matrix: which assays exist per tissue
# ══════════════════════════════════════════════════════════════════════════════
cat("\nBuilding visual table...\n")

# Helper to make a labelled heatmap column strip
make_heat_col <- function(dat, x_col, y_col = "tissue", fill_col, label_col = NULL,
                          fill_lim = NULL, palette = "RdBu",
                          col_title = "", text_size = 2.5) {
  dat2 <- dat %>% mutate(x = !!sym(x_col), y = factor(!!sym(y_col), levels = rev(TISSUE_ORDER)))
  if (is.null(fill_lim)) fill_lim <- max(abs(dat2[[fill_col]]), na.rm = TRUE)
  p <- ggplot(dat2, aes(x = x, y = y, fill = !!sym(fill_col))) +
    geom_tile(colour = "white", linewidth = 0.4)
  if (!is.null(label_col)) {
    p <- p + geom_text(aes(label = !!sym(label_col)), size = text_size, colour = "grey10")
  }
  p + scale_fill_distiller(palette = palette, limits = c(-fill_lim, fill_lim),
                            direction = -1, na.value = "grey88", guide = "none") +
    scale_x_discrete(name = NULL) +
    scale_y_discrete(name = NULL) +
    labs(title = col_title) +
    theme_pub +
    theme(axis.text.x  = element_text(size = 7.5, face = "bold"),
          axis.text.y  = element_text(size = 8),
          plot.title   = element_text(size = 8, face = "bold", hjust = 0.5),
          panel.border = element_blank())
}

# Prepare long-form data for heatmap panels
heat_input <- sex_table %>%
  mutate(tissue = factor(tissue, levels = TISSUE_ORDER),
         lfc_label = sprintf("%.2f", ets_mean))

ets_long <- sex_table %>%
  select(tissue, sex_label, ets_mean) %>%
  mutate(tissue = factor(tissue, levels = TISSUE_ORDER),
         lbl = sprintf("%.2f", ets_mean))

buf_long <- sex_table %>%
  select(tissue, sex_label, buf_mean) %>%
  mutate(tissue = factor(tissue, levels = TISSUE_ORDER),
         lbl = sprintf("%.2f", buf_mean))

ratio_long <- sex_table %>%
  select(tissue, sex_label, buf_ets_ratio) %>%
  mutate(tissue = factor(tissue, levels = TISSUE_ORDER),
         lbl = if_else(is.na(buf_ets_ratio), "—", sprintf("%.2f", buf_ets_ratio)),
         ratio_plot = pmin(pmax(buf_ets_ratio, -3), 3))   # clamp for colour scale

ets_lim <- max(abs(ets_long$ets_mean), na.rm = TRUE)
buf_lim <- max(abs(buf_long$buf_mean), na.rm = TRUE)

# Panel A — ETS heatmap
pA <- ggplot(ets_long, aes(x = sex_label, y = tissue, fill = ets_mean)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = lbl), size = 2.3, colour = "grey10") +
  scale_fill_distiller(palette = "RdBu", limits = c(-ets_lim, ets_lim),
                       direction = -1, na.value = "grey88", guide = "none") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(TISSUE_ORDER)) +
  labs(x = NULL, y = NULL, title = "ETS\nmean log₂FC") +
  theme_pub +
  theme(axis.text.x  = element_text(size = 7.5, face = "bold", colour = c("#2166AC","#D6604D")),
        axis.text.y  = element_text(size = 8.5, face = "plain"),
        plot.title   = element_text(size = 8, face = "bold", hjust = 0.5),
        panel.border = element_blank())

# Panel B — Buffering heatmap
pB <- ggplot(buf_long, aes(x = sex_label, y = tissue, fill = buf_mean)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = lbl), size = 2.3, colour = "grey10") +
  scale_fill_distiller(palette = "RdBu", limits = c(-buf_lim, buf_lim),
                       direction = -1, na.value = "grey88", guide = "none") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(TISSUE_ORDER)) +
  labs(x = NULL, y = NULL, title = "Buffering\nmean log₂FC") +
  theme_pub +
  theme(axis.text.x  = element_text(size = 7.5, face = "bold", colour = c("#2166AC","#D6604D")),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title   = element_text(size = 8, face = "bold", hjust = 0.5),
        panel.border = element_blank())

# Panel C — Buffering/ETS ratio
pC <- ggplot(ratio_long, aes(x = sex_label, y = tissue, fill = ratio_plot)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = lbl), size = 2.2, colour = "grey10") +
  scale_fill_distiller(palette = "PRGn", limits = c(-3, 3),
                       direction = 1, na.value = "grey88", guide = "none") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(TISSUE_ORDER)) +
  labs(x = NULL, y = NULL, title = "Buf/ETS\nratio") +
  theme_pub +
  theme(axis.text.x  = element_text(size = 7.5, face = "bold", colour = c("#2166AC","#D6604D")),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title   = element_text(size = 8, face = "bold", hjust = 0.5),
        panel.border = element_blank())

# Panel D — ETS n_sig counts
nsig_long <- sex_table %>%
  select(tissue, sex_label, ets_n_sig, buf_n_sig) %>%
  pivot_longer(c(ets_n_sig, buf_n_sig), names_to = "metric", values_to = "n_sig") %>%
  mutate(
    tissue = factor(tissue, levels = TISSUE_ORDER),
    col_label = paste0(str_to_title(str_remove(metric,"_n_sig")), "\n", sex_label)
  )

pD <- ggplot(nsig_long, aes(x = col_label, y = tissue, fill = n_sig)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = n_sig), size = 2.3, colour = "grey10") +
  scale_fill_gradient(low = "white", high = "#2C7BB6", na.value = "grey88",
                      guide = "none", limits = c(0, NA)) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(TISSUE_ORDER)) +
  labs(x = NULL, y = NULL, title = "Sig genes\n(adj_p<0.05)") +
  theme_pub +
  theme(axis.text.x  = element_text(size = 6.5, face = "plain"),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title   = element_text(size = 8, face = "bold", hjust = 0.5),
        panel.border = element_blank())

# Panel E — Top Ring 2/3 category + concordance (text columns)
annot_dat <- master_wide %>%
  mutate(
    tissue = factor(tissue, levels = TISSUE_ORDER),
    r23_label    = str_wrap(coalesce(top_r23_cat, "—"), 10),
    gene_label   = if_else(
      !is.na(top_gene_adj_p) & top_gene_adj_p < 0.05,
      paste0(top_gene, "*\n(", top_gene_sex, ")"),
      paste0(coalesce(top_gene,"—"), "\n(", coalesce(top_gene_sex,""), ")")),
    concord_label = if_else(is.na(trx_prot_r_ets), "no\nPROT",
                             sprintf("r=%.2f\n(%s)", trx_prot_r_ets,
                                     if_else(trx_prot_ets_sig > 0, "p<.05", "ns"))),
    coverage_label = assay_coverage
  )

make_text_col <- function(dat, label_col, title, col_fill = "grey95",
                           y_lab = FALSE, text_sz = 2.2) {
  ggplot(dat, aes(x = 1, y = tissue)) +
    geom_tile(fill = col_fill, colour = "white", linewidth = 0.4) +
    geom_text(aes(label = !!sym(label_col)), size = text_sz,
              colour = "grey15", lineheight = 0.9) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_discrete(limits = rev(TISSUE_ORDER)) +
    labs(x = NULL, y = NULL, title = title) +
    theme_pub +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = if (y_lab) element_text(size=8.5) else element_blank(),
          axis.ticks.y = element_blank(),
          plot.title   = element_text(size=8, face="bold", hjust=0.5),
          panel.border = element_blank())
}

pE1 <- make_text_col(annot_dat, "r23_label",     "Top Ring\n2/3 cat.", text_sz = 2.1)
pE2 <- make_text_col(annot_dat, "gene_label",    "Top sig.\ngene*",    text_sz = 2.1)
pE3 <- make_text_col(annot_dat, "concord_label", "Trx-Prot\nconcord.", text_sz = 2.1)
pE4 <- make_text_col(annot_dat, "coverage_label","Assay\ncoverage",    text_sz = 2.4)

# Assemble full figure
fig_table <- (pA | pB | pC | pD | pE1 | pE2 | pE3 | pE4) +
  plot_layout(widths = c(1.8, 1.8, 1.8, 2.4, 1.4, 1.6, 1.4, 1.1)) +
  plot_annotation(
    title    = "Supplementary Table 1 — Tissue Profiles: ETS Expansion vs. Redox Buffering",
    subtitle = paste0(
      "TRNSCRPT, 8-week training vs. sedentary. ",
      "ETS = CI+III+IV (21 genes). Buffering = GSH+TRX+NNT+SOD (15 genes). ",
      "Coverage: P=Protein, Ph=Phospho, A=Acetyl, U=Ubiq, M=Metabolomics. ",
      "* adj_p < 0.05."
    ),
    theme = theme(plot.title    = element_text(size = 10, face = "bold"),
                  plot.subtitle = element_text(size = 7.5, colour = "grey40"))
  )

ggsave("figures/fig_supp_table1_tissue_profiles.pdf",
       fig_table, width = 17, height = 11, device = cairo_pdf)
ggsave("figures/fig_supp_table1_tissue_profiles.png",
       fig_table, width = 17, height = 11, dpi = 200)
cat("Saved figures/fig_supp_table1_tissue_profiles (.pdf/.png)\n")

# ══════════════════════════════════════════════════════════════════════════════
# BLOCK 8 — Sex-averaged compact summary (simpler CSV for quick reference)
# ══════════════════════════════════════════════════════════════════════════════
avg_summary <- master_csv %>%
  mutate(
    # Interpret buffering/ETS ratio narrative
    redox_status = case_when(
      ets_mean_avg > 0.10 & buf_mean_avg > 0.05 &
        (buf_ets_ratio_male > 0.5 | buf_ets_ratio_female > 0.5) ~ "Balanced expansion",
      ets_mean_avg > 0.10 & buf_mean_avg < 0.05                  ~ "ETS-only (deficit)",
      ets_mean_avg < -0.10                                        ~ "ETS suppression",
      abs(ets_mean_avg) <= 0.10 & abs(buf_mean_avg) <= 0.05      ~ "Minimal response",
      TRUE                                                         ~ "Mixed / moderate"
    ),
    # Sex difference flag
    sex_diff_ets = abs(ets_mean_male - ets_mean_female) > 0.15
  ) %>%
  select(tissue, ets_mean_avg, buf_mean_avg, buf_ets_ratio_male, buf_ets_ratio_female,
         redox_status, sex_diff_ets, top_r23_cat, top_gene, top_gene_adj_p,
         trx_prot_r_ets, assay_coverage)

write_csv(avg_summary, "data/tissue_profiles_master_avg.csv")
cat("Saved data/tissue_profiles_master_avg.csv\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# Print full console summary
# ══════════════════════════════════════════════════════════════════════════════
cat("=== SUPPLEMENTARY TABLE 1 — FULL SUMMARY ===\n\n")
avg_summary %>%
  arrange(match(tissue, TISSUE_ORDER)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\n=== REDOX STATUS DISTRIBUTION ===\n")
print(table(avg_summary$redox_status))

cat("\n=== TISSUES WITH LARGEST ETS/BUFFERING GAP ===\n")
avg_summary %>%
  filter(!is.na(buf_ets_ratio_male)) %>%
  arrange(buf_ets_ratio_male) %>%
  select(tissue, ets_mean_avg, buf_mean_avg, buf_ets_ratio_male, redox_status) %>%
  head(10) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\n=== TISSUES WITH SEX DIFFERENCES IN ETS EXPANSION ===\n")
master_csv %>%
  filter(abs(ets_mean_male - ets_mean_female) > 0.15) %>%
  select(tissue, ets_mean_male, ets_mean_female, ets_n_sig_male, ets_n_sig_female) %>%
  arrange(desc(abs(ets_mean_male - ets_mean_female))) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\nDone.\n")
