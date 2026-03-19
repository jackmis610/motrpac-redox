# =============================================================================
# 03_figures.R
# Purpose : Generate primary display figures for the Neufer framework analysis.
#           Fig 1 — Tissue × functional-group heatmap (8w TRNSCRPT)
#           Fig 2 — BAT deep dive: ETS gene trajectories + metabolomics
# Inputs  : data/neufer_da_all.csv   (from 02_query_neufer.R)
#           MotrpacRatTraining6moData (METAB_BAT_DA for metabolomics panel)
# Outputs : figures/fig1_tissue_category_heatmap.png/.pdf
#           figures/fig2a_bat_gene_heatmap.png/.pdf
#           figures/fig2bcd_bat_trajectories_metabolomics.png/.pdf
#           data/bat_gsh_discordance.csv
# =============================================================================
# 03: Neufer framework figures
#   Fig 1 — Tissue × functional-group heatmap (8w, TRNSCRPT)
#   Fig 2 — BAT deep dive: gene trajectories + metabolomics
# ============================================================

suppressPackageStartupMessages({
  library(MotrpacRatTraining6moData)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

dir.create("figures", showWarnings = FALSE)

da <- read_csv("data/neufer_da_all.csv", show_col_types = FALSE)

# ── Shared aesthetics ──────────────────────────────────────────────────────────
DIVPAL  <- "RdBu"       # diverging: blue = down, red = up
SEX_COL <- c(Male = "#2166AC", Female = "#D6604D")

# Framework category display order (Neufer circuit logic: ETS → PMF → NADPH → buffers → fuels)
CAT_ORDER <- c("ets_CI","ets_CIII","ets_CIV","cyt_c",
               "atp_synthase","nnt","nadph_other",
               "gsh","trx","sod",
               "q_pool","beta_ox","pdh","shuttles","ucp")

CAT_LABELS <- c(
  ets_CI       = "Complex I",
  ets_CIII     = "Complex III",
  ets_CIV      = "Complex IV",
  cyt_c        = "Cyt c",
  atp_synthase = "ATP Synthase",
  nnt          = "NNT",
  nadph_other  = "NADPH (other)",
  gsh          = "GSH system",
  trx          = "TRX system",
  sod          = "SOD",
  q_pool       = "Q-pool feeders",
  beta_ox      = "β-oxidation",
  pdh          = "PDH complex",
  shuttles     = "Shuttles",
  ucp          = "UCPs"
)

# Tissue display order (skeletal muscle first, then oxidative tissues, then others)
TISSUE_ORDER <- c("SKM-GN","SKM-VL","HEART","BAT","WAT-SC",
                  "LIVER","KIDNEY","LUNG","ADRNL",
                  "BLOOD","CORTEX","HIPPOC","HYPOTH","VENACV",
                  "COLON","SMLINT","SPLEEN","OVARY","TESTES")

theme_clean <- theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text       = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank(),
    legend.key.size  = unit(0.4, "cm")
  )

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Tissue × functional-group heatmap (TRNSCRPT, 8w)
# ══════════════════════════════════════════════════════════════════════════════
cat("Building Figure 1: tissue × category heatmap...\n")

trx8 <- da %>%
  filter(assay == "TRNSCRPT", comparison_group == "8w") %>%
  mutate(
    framework_category = factor(framework_category, levels = CAT_ORDER),
    tissue = factor(tissue, levels = TISSUE_ORDER)
  )

# Mean logFC per tissue × category × sex; flag if ≥1 gene adj_p < 0.05
heat_dat <- trx8 %>%
  group_by(tissue, framework_category, sex) %>%
  summarise(
    mean_logFC  = mean(logFC, na.rm = TRUE),
    n_genes     = n_distinct(gene_symbol),
    n_sig       = sum(adj_p_value < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sig_label = case_when(
      n_sig >= 3 ~ "***",
      n_sig == 2 ~ "**",
      n_sig == 1 ~ "*",
      TRUE       ~ ""
    ),
    sex_label = str_to_title(sex)
  )

lim <- max(abs(heat_dat$mean_logFC), na.rm = TRUE)
lim <- ceiling(lim * 10) / 10   # round up to nearest 0.1

fig1 <- ggplot(heat_dat, aes(x = framework_category, y = tissue, fill = mean_logFC)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sig_label), size = 2.8, vjust = 0.75, colour = "grey20") +
  scale_fill_distiller(
    palette  = DIVPAL,
    limits   = c(-lim, lim),
    direction = -1,
    name     = "Mean log₂FC\nvs. sedentary",
    na.value = "grey85"
  ) +
  scale_x_discrete(labels = CAT_LABELS) +
  scale_y_discrete(limits = rev(levels(heat_dat$tissue))) +
  facet_wrap(~ sex_label, ncol = 1) +
  labs(
    title    = "Neufer Framework Gene Expression — 8-Week Training Response",
    subtitle = "Mean log₂FC (TRNSCRPT) vs. sedentary control; * ≥1 gene adj_p < 0.05",
    x        = NULL,
    y        = NULL
  ) +
  theme_clean +
  theme(
    axis.text.x      = element_text(angle = 40, hjust = 1, size = 8.5),
    axis.text.y      = element_text(size = 8.5),
    panel.border     = element_blank(),
    legend.position  = "right"
  )

ggsave("figures/fig1_tissue_category_heatmap.pdf",
       fig1, width = 11, height = 10, device = cairo_pdf)
ggsave("figures/fig1_tissue_category_heatmap.png",
       fig1, width = 11, height = 10, dpi = 200)
cat("  Saved fig1_tissue_category_heatmap (.pdf/.png)\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — BAT deep dive
# Panel A: Gene heatmap (TRNSCRPT, all timepoints)
# Panel B: Category-level temporal trajectories
# Panel C: Redox metabolomics (METAB_BAT_DA)
# ══════════════════════════════════════════════════════════════════════════════
cat("\nBuilding Figure 2: BAT deep dive...\n")

# ── 2A: Gene-level heatmap ────────────────────────────────────────────────────
bat_trx <- da %>%
  filter(tissue == "BAT", assay == "TRNSCRPT") %>%
  mutate(
    framework_category = factor(framework_category, levels = CAT_ORDER),
    comparison_group   = factor(comparison_group,   levels = c("1w","2w","4w","8w")),
    sex_label          = str_to_title(sex)
  )

# Flag significance per gene × timepoint × sex
bat_trx <- bat_trx %>%
  mutate(sig_dot = if_else(adj_p_value < 0.05, "●", ""))

# Gene ordering: by category (as in CAT_ORDER) then alphabetical within category
gene_order <- bat_trx %>%
  distinct(gene_symbol, framework_category) %>%
  arrange(framework_category, gene_symbol) %>%
  pull(gene_symbol)

bat_trx <- bat_trx %>%
  mutate(gene_symbol = factor(gene_symbol, levels = gene_order))

# Category colour bar (for left annotation)
cat_cols <- c(
  ets_CI       = "#1F78B4", ets_CIII  = "#33A02C", ets_CIV   = "#B2DF8A",
  cyt_c        = "#6A3D9A", atp_synthase = "#CAB2D6", nnt     = "#E31A1C",
  nadph_other  = "#FF7F00", gsh       = "#FDBF6F",   trx     = "#A6CEE3",
  sod          = "#B15928", q_pool    = "#FB9A99",   beta_ox  = "#FFFF99",
  pdh          = "#1A1A1A", shuttles  = "#999999",   ucp     = "#D9D9D9"
)

lim_bat <- max(abs(bat_trx$logFC), na.rm = TRUE)
lim_bat <- ceiling(lim_bat * 5) / 5   # round to nearest 0.2

fig2a <- ggplot(bat_trx, aes(x = comparison_group, y = gene_symbol, fill = logFC)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  geom_text(aes(label = sig_dot), size = 2.5, colour = "black", vjust = 0.75) +
  scale_fill_distiller(
    palette   = DIVPAL,
    limits    = c(-lim_bat, lim_bat),
    direction = -1,
    name      = "log₂FC",
    na.value  = "grey85"
  ) +
  scale_y_discrete(limits = rev(gene_order)) +
  facet_grid(framework_category ~ sex_label,
             scales = "free_y", space = "free_y",
             labeller = labeller(framework_category = CAT_LABELS)) +
  labs(
    title = "BAT — Neufer Gene Temporal Response (TRNSCRPT)",
    subtitle = "● adj_p < 0.05",
    x = "Training duration",
    y = NULL
  ) +
  theme_clean +
  theme(
    axis.text.y      = element_text(size = 6.5, face = "italic"),
    axis.text.x      = element_text(size = 8),
    strip.text.y     = element_text(size = 6.5, angle = 0, hjust = 0),
    strip.text.x     = element_text(size = 9),
    panel.border     = element_rect(colour = "grey80", fill = NA),
    legend.position  = "right"
  )

# ── 2B: Category-level temporal trajectories ─────────────────────────────────
cat_traj <- bat_trx %>%
  group_by(framework_category, comparison_group, sex_label) %>%
  summarise(
    mean_lfc = mean(logFC, na.rm = TRUE),
    se_lfc   = sd(logFC, na.rm = TRUE) / sqrt(n_distinct(gene_symbol)),
    .groups  = "drop"
  )

fig2b <- ggplot(cat_traj,
                aes(x    = comparison_group,
                    y    = mean_lfc,
                    colour = framework_category,
                    group  = framework_category)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_ribbon(aes(ymin  = mean_lfc - se_lfc,
                  ymax  = mean_lfc + se_lfc,
                  fill  = framework_category),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_colour_manual(values = cat_cols, labels = CAT_LABELS, name = "Category") +
  scale_fill_manual(  values = cat_cols, labels = CAT_LABELS, name = "Category") +
  facet_wrap(~ sex_label, ncol = 2) +
  labs(
    title    = "BAT — Category Mean log₂FC Trajectories",
    subtitle = "Mean ± SE across genes within category (TRNSCRPT)",
    x        = "Training duration",
    y        = "Mean log₂FC"
  ) +
  theme_clean +
  theme(legend.text = element_text(size = 7.5),
        legend.key.size = unit(0.35, "cm"))

# ── 2C: Redox metabolomics ────────────────────────────────────────────────────
data(METAB_BAT_DA)

# Key metabolites for the Neufer story, grouped by theme
metab_groups <- tribble(
  ~feature_ID,             ~metab_group,
  "NAD+",                  "NAD redox",
  "NADH",                  "NAD redox",
  "NADP+",                 "NADP redox",
  "NADPH",                 "NADP redox",
  "FAD",                   "Flavin",
  "glutathione",           "Glutathione",
  "glutathione reduced",   "Glutathione",
  "oxidized glutathione",  "Glutathione",
  "succinate",             "TCA intermediates",
  "fumarate",              "TCA intermediates",
  "malate",                "TCA intermediates",
  "pyruvate",              "TCA intermediates",
  "lactate",               "TCA intermediates"
)

metab_order <- c("NAD redox","NADP redox","Flavin","Glutathione","TCA intermediates")

bat_metab <- METAB_BAT_DA %>%
  inner_join(metab_groups, by = "feature_ID") %>%
  mutate(
    comparison_group = factor(comparison_group, levels = c("1w","2w","4w","8w")),
    metab_group      = factor(metab_group, levels = metab_order),
    sex_label        = str_to_title(sex),
    sig_dot          = if_else(adj_p_value < 0.05, "●", "")
  ) %>%
  # Average duplicate feature_IDs (some metabolites have 2 platform entries)
  group_by(feature_ID, metab_group, comparison_group, sex_label, sig_dot) %>%
  summarise(logFC = mean(logFC, na.rm = TRUE),
            adj_p_value = min(adj_p_value, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(sig_dot = if_else(adj_p_value < 0.05, "●", ""))

fig2c <- ggplot(bat_metab,
                aes(x = comparison_group, y = logFC,
                    group = interaction(feature_ID, sex_label),
                    colour = sex_label, linetype = feature_ID)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 0.75) +
  geom_point(aes(shape = sex_label), size = 2.2) +
  geom_text(aes(label = sig_dot), nudge_x = 0.15, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = SEX_COL, name = "Sex") +
  scale_shape_manual( values = c(Male = 16, Female = 17), name = "Sex") +
  facet_wrap(~ metab_group, scales = "free_y", ncol = 2) +
  labs(
    title    = "BAT — Redox Metabolomics (METAB_BAT_DA)",
    subtitle = "● adj_p < 0.05; log₂FC vs. sedentary",
    x        = "Training duration",
    y        = "log₂FC"
  ) +
  theme_clean +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(size = 8)
  ) +
  guides(linetype = guide_legend(title = "Metabolite", ncol = 2, override.aes = list(size = 0.8)))

# ── 2D: Multi-omic discordance — GSH system transcript vs. metabolite ─────────
# Flag: genes where TRNSCRPT logFC > 0.2 at 8w but corresponding metabolite ≤ 0

gsh_trx <- bat_trx %>%
  filter(framework_category == "gsh", comparison_group == "8w") %>%
  select(gene_symbol, sex_label, trx_lfc = logFC, trx_adj_p = adj_p_value)

# Glutathione metabolites at 8w
gsh_metab_8w <- bat_metab %>%
  filter(metab_group == "Glutathione", comparison_group == "8w") %>%
  group_by(sex_label) %>%
  summarise(
    mean_gsh_lfc = mean(logFC[feature_ID %in% c("glutathione","glutathione reduced")], na.rm = TRUE),
    mean_gssg_lfc = mean(logFC[feature_ID == "oxidized glutathione"], na.rm = TRUE),
    .groups = "drop"
  )

discordance <- gsh_trx %>%
  left_join(gsh_metab_8w, by = "sex_label") %>%
  mutate(
    discordant = case_when(
      trx_lfc > 0.15 & mean_gsh_lfc < 0  ~ "transcript ↑ / GSH ↓",
      trx_lfc < -0.15 & mean_gsh_lfc > 0 ~ "transcript ↓ / GSH ↑",
      TRUE                                  ~ "concordant"
    )
  )

fig2d <- ggplot(discordance,
                aes(x = trx_lfc, y = mean_gsh_lfc,
                    colour = discordant, label = gene_symbol)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_point(size = 3, alpha = 0.85) +
  ggrepel::geom_text_repel(size = 2.8, show.legend = FALSE,
                            max.overlaps = 15, seed = 42) +
  scale_colour_manual(
    values = c("concordant" = "grey50",
               "transcript ↑ / GSH ↓" = "#D6604D",
               "transcript ↓ / GSH ↑" = "#4393C3"),
    name   = "Discordance"
  ) +
  facet_wrap(~ sex_label) +
  labs(
    title    = "BAT — GSH System: 8w Transcript vs. Metabolite (Multi-omic Discordance)",
    subtitle = "x-axis: gene log₂FC (TRNSCRPT 8w); y-axis: mean GSH pool log₂FC (metabolomics 8w)",
    x        = "Gene log₂FC (TRNSCRPT, 8w)",
    y        = "Mean GSH pool log₂FC (metabolomics, 8w)"
  ) +
  theme_clean +
  theme(legend.position = "bottom")

# ── Assemble and save Figure 2 ────────────────────────────────────────────────

# Check if ggrepel is available; if not fall back
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  warning("ggrepel not found — discordance labels will use geom_text instead")
  fig2d <- fig2d +
    geom_text(size = 2.5, vjust = -0.5, show.legend = FALSE)
}

# Save panels individually (gene heatmap is tall; compose separately)
ggsave("figures/fig2a_bat_gene_heatmap.pdf",
       fig2a, width = 9, height = 18, device = cairo_pdf)
ggsave("figures/fig2a_bat_gene_heatmap.png",
       fig2a, width = 9, height = 18, dpi = 180)
cat("  Saved fig2a_bat_gene_heatmap (.pdf/.png)\n")

fig2bcd <- (fig2b / fig2c / fig2d) +
  plot_annotation(
    title = "Figure 2 — BAT Deep Dive: Temporal Trajectories & Redox Metabolomics",
    tag_levels = list(c("B","C","D"))
  )

ggsave("figures/fig2bcd_bat_trajectories_metabolomics.pdf",
       fig2bcd, width = 12, height = 18, device = cairo_pdf)
ggsave("figures/fig2bcd_bat_trajectories_metabolomics.png",
       fig2bcd, width = 12, height = 18, dpi = 180)
cat("  Saved fig2bcd_bat_trajectories_metabolomics (.pdf/.png)\n")

# ── Discordance summary table ─────────────────────────────────────────────────
cat("\n=== GSH multi-omic discordance at 8w (BAT) ===\n")
disc_summary <- discordance %>%
  arrange(sex_label, desc(abs(trx_lfc))) %>%
  select(sex_label, gene_symbol, trx_lfc, trx_adj_p, mean_gsh_lfc, mean_gssg_lfc, discordant)
print(disc_summary, n = 30)

write_csv(disc_summary, "data/bat_gsh_discordance.csv")
cat("\nSaved data/bat_gsh_discordance.csv\n")

# ── Summary stats ─────────────────────────────────────────────────────────────
cat("\n=== Figure 1 summary ===\n")
top_tiles <- heat_dat %>%
  filter(!is.na(mean_logFC)) %>%
  slice_max(abs(mean_logFC), n = 10) %>%
  select(tissue, framework_category, sex, mean_logFC, n_sig)
print(top_tiles)

cat("\n=== BAT significant genes at any timepoint ===\n")
bat_sig <- bat_trx %>%
  filter(adj_p_value < 0.05) %>%
  distinct(gene_symbol, framework_category, sex_label, comparison_group) %>%
  arrange(framework_category, gene_symbol)
print(bat_sig, n = 40)

cat("\nDone — all figures saved to figures/\n")
