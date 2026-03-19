# =============================================================================
# 04_sensitivity.R
# Purpose : Statistical rigor analyses for the ETS–buffering scaling claim:
#           (1) GPx4 tissue specificity across the ETS expansion gradient
#           (2) Transcript–protein concordance for ETS subunits and NNT
#           (3) ETS vs. buffering regression (the main paper figure, Fig 2)
#           Cook's distance and leave-one-out diagnostics are in
#           supp_cooks_distance.R; BAT-excluded sensitivity in
#           supp_sensitivity_noBat.R.
# Inputs  : data/neufer_da_all.csv   (from 02_query_neufer.R)
# Outputs : figures/fig3a_gpx4_tissue_specificity.png/.pdf
#           figures/fig3b_trx_prot_concordance.png/.pdf
#           figures/fig3c_ets_vs_buffering_correlation.png/.pdf
#           figures/fig2_ets_buffering_v4.png/.pdf  (publication-ready version)
#           data/gpx4_tissue_specificity.csv
#           data/trx_prot_concordance.csv
#           data/ets_buffering_correlation_data.csv
#           data/ets_buffering_regression_stats.csv
#           data/regression_sensitivity_table.csv
# =============================================================================
# 04: Statistical rigor — three analyses
#   (1) GPx4 tissue specificity across ETS expansion gradient
#   (2) Transcript-protein concordance for ETS subunits + NNT
#   (3) ETS vs. buffering circuit correlation — the paper figure
# ============================================================

suppressPackageStartupMessages({
  library(MotrpacRatTraining6moData)
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
})

dir.create("figures", showWarnings = FALSE)
dir.create("data",    showWarnings = FALSE)

da <- read_csv("data/neufer_da_all.csv", show_col_types = FALSE)

SEX_COL   <- c(Male = "#2166AC", Female = "#D6604D")
SEX_SHAPE <- c(Male = 16L,       Female = 17L)

theme_pub <- theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text       = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank(),
    legend.key.size  = unit(0.4, "cm")
  )

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 1 — GPx4 tissue specificity
# Does Gpx4 downregulation track with ETS expansion, or is BAT unique?
# ══════════════════════════════════════════════════════════════════════════════
cat("=== Analysis 1: GPx4 tissue specificity ===\n\n")

# ETS expansion per tissue × sex at 8w (TRNSCRPT)
ets_8w <- da %>%
  filter(assay == "TRNSCRPT", comparison_group == "8w",
         framework_category %in% c("ets_CI","ets_CIII","ets_CIV")) %>%
  group_by(tissue, sex) %>%
  summarise(
    ets_mean_lfc = mean(logFC, na.rm = TRUE),
    ets_n_sig    = sum(adj_p_value < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

# GPx4 TRNSCRPT at 8w per tissue × sex
gpx4_trx <- da %>%
  filter(assay == "TRNSCRPT", comparison_group == "8w", gene_symbol == "Gpx4") %>%
  select(tissue, sex, gpx4_lfc = logFC, gpx4_adj_p = adj_p_value, gpx4_p = p_value)

# GPx4 PROT at 8w (where available)
gpx4_prot <- da %>%
  filter(assay == "PROT", comparison_group == "8w", gene_symbol == "Gpx4") %>%
  select(tissue, sex, gpx4_prot_lfc = logFC, gpx4_prot_adj_p = adj_p_value)

# Join
gpx4_dat <- ets_8w %>%
  left_join(gpx4_trx,  by = c("tissue","sex")) %>%
  left_join(gpx4_prot, by = c("tissue","sex")) %>%
  mutate(
    sex_label = str_to_title(sex),
    trx_sig   = case_when(
      gpx4_adj_p < 0.001 ~ "***",
      gpx4_adj_p < 0.01  ~ "**",
      gpx4_adj_p < 0.05  ~ "*",
      TRUE                ~ ""
    ),
    point_size = if_else(trx_sig != "", 3.5, 2.5)
  )

# Correlation across tissues (within sex)
cor_stats <- gpx4_dat %>%
  group_by(sex_label) %>%
  summarise(
    r    = cor(ets_mean_lfc, gpx4_lfc, use = "complete.obs"),
    n    = sum(!is.na(gpx4_lfc)),
    pval = cor.test(ets_mean_lfc, gpx4_lfc)$p.value,
    .groups = "drop"
  )
cat("GPx4 ~ ETS expansion correlation (TRNSCRPT, 8w):\n")
print(cor_stats)

# Print GPx4 summary across tissues
cat("\nGPx4 TRNSCRPT 8w — significant hits:\n")
gpx4_dat %>%
  filter(trx_sig != "") %>%
  select(tissue, sex_label, ets_mean_lfc, gpx4_lfc, gpx4_adj_p, gpx4_prot_lfc) %>%
  arrange(gpx4_adj_p) %>%
  print()

cat("\nGPx4 TRNSCRPT 8w — all tissues ranked by ETS expansion:\n")
gpx4_dat %>%
  arrange(sex, desc(ets_mean_lfc)) %>%
  select(tissue, sex_label, ets_mean_lfc, ets_n_sig, gpx4_lfc, gpx4_adj_p) %>%
  print(n=40)

# Fig 3A — scatter: ETS expansion vs GPx4 logFC
fig3a <- ggplot(gpx4_dat,
                aes(x = ets_mean_lfc, y = gpx4_lfc,
                    colour = sex_label, shape = sex_label,
                    label  = tissue)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.15) +
  geom_point(aes(size = trx_sig != ""), show.legend = TRUE) +
  geom_text_repel(size = 2.8, max.overlaps = 20, seed = 42,
                  segment.size = 0.3, min.segment.length = 0.2) +
  # annotate significant points
  geom_text(aes(label = trx_sig), vjust = -1.0, size = 3.5,
            show.legend = FALSE, colour = "grey20") +
  geom_text(
    data = cor_stats,
    aes(x = -Inf, y = Inf,
        label = sprintf("r = %.2f, p = %.3f", r, pval),
        colour = sex_label),
    hjust = -0.1, vjust = 1.5, size = 3, inherit.aes = FALSE
  ) +
  scale_colour_manual(values = SEX_COL, name = "Sex") +
  scale_shape_manual( values = SEX_SHAPE, name = "Sex") +
  scale_size_manual(  values = c(`FALSE` = 2.2, `TRUE` = 3.8),
                      guide  = "none") +
  facet_wrap(~ sex_label) +
  labs(
    title    = "GPx4 Expression vs. ETS Expansion Across Tissues (TRNSCRPT, 8w)",
    subtitle = "* adj_p < 0.05  ** adj_p < 0.01  *** adj_p < 0.001",
    x        = "Mean ETS log\u2082FC (CI+III+IV genes)",
    y        = "Gpx4 log\u2082FC"
  ) +
  theme_pub +
  theme(legend.position = "none")

# Fig 3A-inset: GPx4 PROT vs TRNSCRPT for tissues with protein data
gpx4_prot_cmp <- gpx4_dat %>%
  filter(!is.na(gpx4_prot_lfc)) %>%
  mutate(sex_label = str_to_title(sex))

fig3a_inset <- ggplot(gpx4_prot_cmp,
                      aes(x = gpx4_lfc, y = gpx4_prot_lfc,
                          colour = sex_label, label = tissue)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_point(size = 2.5) +
  geom_text_repel(size = 2.5, seed = 42, max.overlaps = 15,
                  segment.size = 0.25) +
  scale_colour_manual(values = SEX_COL, name = "Sex") +
  labs(
    title    = "GPx4: Transcript vs. Protein (8w)",
    subtitle = "Dashed = identity (1:1); tissues with PROT data only",
    x        = "Gpx4 log\u2082FC (TRNSCRPT)",
    y        = "Gpx4 log\u2082FC (PROT)"
  ) +
  theme_pub +
  theme(legend.position = "right")

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 2 — Transcript-protein concordance (7 tissues)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n\n=== Analysis 2: Transcript-protein concordance ===\n\n")

PROT_TISSUES <- c("CORTEX","HEART","KIDNEY","LIVER","LUNG","SKM-GN","WAT-SC")

trx_8w <- da %>%
  filter(assay == "TRNSCRPT", comparison_group == "8w",
         tissue %in% PROT_TISSUES) %>%
  select(gene_symbol, framework_category, tissue, sex, trx_lfc = logFC,
         trx_adj_p = adj_p_value)

prot_8w <- da %>%
  filter(assay == "PROT", comparison_group == "8w") %>%
  select(gene_symbol, tissue, sex, prot_lfc = logFC, prot_adj_p = adj_p_value)

trx_prot <- trx_8w %>%
  inner_join(prot_8w, by = c("gene_symbol","tissue","sex"))

# Per-tissue Pearson r (all Neufer genes)
concordance <- trx_prot %>%
  group_by(tissue, sex) %>%
  summarise(
    r_all  = cor(trx_lfc, prot_lfc, use = "complete.obs"),
    p_all  = cor.test(trx_lfc, prot_lfc)$p.value,
    n_all  = n(),
    .groups = "drop"
  )

# ETS-specific concordance
concordance_ets <- trx_prot %>%
  filter(framework_category %in% c("ets_CI","ets_CIII","ets_CIV")) %>%
  group_by(tissue, sex) %>%
  summarise(
    r_ets = cor(trx_lfc, prot_lfc, use = "complete.obs"),
    p_ets = cor.test(trx_lfc, prot_lfc)$p.value,
    n_ets = n(),
    .groups = "drop"
  )

concordance_full <- concordance %>%
  left_join(concordance_ets, by = c("tissue","sex")) %>%
  mutate(sex_label = str_to_title(sex)) %>%
  arrange(sex, desc(r_ets))

cat("Transcript-protein concordance (Pearson r) at 8w:\n")
print(concordance_full %>%
      select(tissue, sex_label, r_all, p_all, r_ets, p_ets, n_ets), n=20)

# NNT protein in tissues with PROT data
cat("\nNNT TRNSCRPT vs PROT at 8w:\n")
trx_prot %>%
  filter(gene_symbol == "Nnt") %>%
  select(tissue, sex, trx_lfc, trx_adj_p, prot_lfc, prot_adj_p) %>%
  arrange(sex, tissue) %>%
  print()

# NNT acetylation (HEART, LIVER) at 8w
cat("\nNNT acetylation (ACETYL) at 8w — all sites:\n")
nnt_acetyl <- da %>%
  filter(gene_symbol == "Nnt", assay == "ACETYL", comparison_group == "8w") %>%
  select(tissue, sex, feature_ID, logFC, adj_p_value, p_value) %>%
  arrange(adj_p_value)
print(nnt_acetyl, n = 30)

# ETS acetylation: any significant sites at 8w?
cat("\nETS acetylation significant sites at 8w (adj_p < 0.05):\n")
ets_acetyl_sig <- da %>%
  filter(framework_category %in% c("ets_CI","ets_CIII","ets_CIV"),
         assay == "ACETYL", comparison_group == "8w",
         adj_p_value < 0.05) %>%
  select(gene_symbol, framework_category, tissue, sex, feature_ID, logFC, adj_p_value) %>%
  arrange(adj_p_value)
print(ets_acetyl_sig, n = 30)

# Fig 3B — scatter: TRNSCRPT vs PROT logFC for ETS genes across 7 tissues
trx_prot_ets <- trx_prot %>%
  filter(framework_category %in% c("ets_CI","ets_CIII","ets_CIV")) %>%
  mutate(sex_label = str_to_title(sex))

fig3b <- ggplot(trx_prot_ets,
                aes(x = trx_lfc, y = prot_lfc,
                    colour = tissue)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey70", linewidth = 0.3) +
  geom_smooth(aes(group = tissue), method = "lm", se = FALSE,
              linewidth = 0.4, alpha = 0.5) +
  geom_point(size = 2.2, alpha = 0.8) +
  # per-panel r annotation
  geom_text(
    data = concordance_ets %>% mutate(sex_label = str_to_title(sex)),
    aes(x = -Inf, y = Inf,
        label = sprintf("r = %.2f", r_ets),
        colour = tissue),
    hjust = -0.1, vjust = 1.8, size = 2.8, inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  scale_colour_brewer(palette = "Dark2", name = "Tissue") +
  facet_wrap(~ sex_label) +
  labs(
    title    = "ETS Subunit Transcript vs. Protein Concordance (8w)",
    subtitle = "Dashed = 1:1 identity; each point = one ETS gene in one tissue",
    x        = "log\u2082FC (TRNSCRPT)",
    y        = "log\u2082FC (PROT)"
  ) +
  theme_pub

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 3 — ETS vs. buffering circuit correlation (the paper figure)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n\n=== Analysis 3: ETS vs. buffering circuit correlation ===\n\n")

# Compute per tissue × sex at 8w, TRNSCRPT
ets_categories     <- c("ets_CI","ets_CIII","ets_CIV","cyt_c","atp_synthase")
buffering_categories <- c("gsh","trx","nnt","sod")

trx_8w_all <- da %>%
  filter(assay == "TRNSCRPT", comparison_group == "8w")

ets_means <- trx_8w_all %>%
  filter(framework_category %in% ets_categories) %>%
  group_by(tissue, sex) %>%
  summarise(
    ets_mean  = mean(logFC, na.rm = TRUE),
    ets_se    = sd(logFC,   na.rm = TRUE) / sqrt(n_distinct(gene_symbol)),
    ets_n_sig = sum(adj_p_value < 0.05, na.rm = TRUE),
    ets_n     = n_distinct(gene_symbol),
    .groups = "drop"
  )

buf_means <- trx_8w_all %>%
  filter(framework_category %in% buffering_categories) %>%
  group_by(tissue, sex) %>%
  summarise(
    buf_mean  = mean(logFC, na.rm = TRUE),
    buf_se    = sd(logFC,   na.rm = TRUE) / sqrt(n_distinct(gene_symbol)),
    buf_n_sig = sum(adj_p_value < 0.05, na.rm = TRUE),
    buf_n     = n_distinct(gene_symbol),
    .groups = "drop"
  )

fig3c_dat <- ets_means %>%
  inner_join(buf_means, by = c("tissue","sex")) %>%
  mutate(sex_label = str_to_title(sex))

# Formal regression per sex: buf ~ ets
# H0: slope = 1 (1:1 scaling); H1: slope < 1 (buffering lags ETS)
reg_stats <- fig3c_dat %>%
  group_by(sex_label) %>%
  summarise({
    fit  <- lm(buf_mean ~ ets_mean)
    cf   <- coef(fit)
    se   <- summary(fit)$coefficients["ets_mean","Std. Error"]
    r2   <- summary(fit)$r.squared
    # one-sided t-test: slope < 1
    t_vs1 <- (cf["ets_mean"] - 1) / se
    p_vs1 <- pt(t_vs1, df = n() - 2)  # one-tailed lower
    tibble(
      intercept = cf["(Intercept)"],
      slope     = cf["ets_mean"],
      slope_se  = se,
      r2        = r2,
      p_slope_lt1 = p_vs1,          # p(slope < 1)
      p_slope_any = 2 * pt(-abs(t_vs1), df = n() - 2),  # two-tailed
      n         = n()
    )
  }, .groups = "drop")

cat("ETS vs. buffering regression (TRNSCRPT, 8w):\n")
print(reg_stats)

cat("\nData points:\n")
fig3c_dat %>%
  arrange(sex, desc(ets_mean)) %>%
  select(tissue, sex_label, ets_mean, buf_mean, ets_n_sig, buf_n_sig) %>%
  print(n = 40)

# Annotation strings for plot
reg_labels <- reg_stats %>%
  mutate(
    lab = sprintf(
      "slope = %.2f (SE %.2f)\nR\u00b2 = %.2f\np(slope < 1) = %.3f",
      slope, slope_se, r2, p_slope_lt1
    )
  )

# Fig 3C — the paper figure
fig3c <- ggplot(fig3c_dat,
                aes(x     = ets_mean,
                    y     = buf_mean,
                    colour = sex_label,
                    shape  = sex_label,
                    label  = tissue)) +
  # Reference: 1:1 line
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey40",
              linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dotted",
             colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted",
             colour = "grey70", linewidth = 0.3) +
  # Regression lines (inherit.aes = FALSE avoids label-in-smooth warning)
  geom_smooth(aes(x = ets_mean, y = buf_mean, colour = sex_label, fill = sex_label),
              method = "lm", se = TRUE, linewidth = 1.0,
              alpha = 0.15, inherit.aes = FALSE) +
  # Error bars (SE across genes within category)
  geom_errorbar(aes(ymin = buf_mean - buf_se, ymax = buf_mean + buf_se),
                width = 0, linewidth = 0.35, alpha = 0.5) +
  geom_errorbar(aes(xmin = ets_mean - ets_se, xmax = ets_mean + ets_se),
                width = 0, linewidth = 0.35, alpha = 0.5,
                orientation = "y") +
  geom_point(size = 3.2) +
  geom_text_repel(size = 2.7, max.overlaps = 20, seed = 42,
                  segment.size = 0.25, min.segment.length = 0.3,
                  show.legend = FALSE) +
  # Regression annotation
  geom_text(
    data = reg_labels,
    aes(x = -Inf, y = Inf, label = lab, colour = sex_label),
    hjust = -0.05, vjust = 1.3, size = 2.8,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  scale_colour_manual(values = SEX_COL,   name = "Sex") +
  scale_shape_manual( values = SEX_SHAPE, name = "Sex") +
  facet_wrap(~ sex_label) +
  labs(
    title    = "ETS Expansion vs. Redox Buffering Capacity Across 19 Tissues",
    subtitle = "TRNSCRPT, 8-week training; dashed = 1:1 (matched scaling); SE bars = gene-level variance",
    x        = "Mean log\u2082FC — ETS + ATP synthase (OXPHOS expansion)",
    y        = "Mean log\u2082FC — GSH + TRX + NNT + SOD (redox buffering)"
  ) +
  theme_pub +
  theme(legend.position = "bottom")

# ══════════════════════════════════════════════════════════════════════════════
# Save figures
# ══════════════════════════════════════════════════════════════════════════════
fig3_top <- (fig3a | fig3a_inset) +
  plot_layout(widths = c(2, 1)) +
  plot_annotation(tag_levels = list(c("A","A\u2019")))

ggsave("figures/fig3a_gpx4_tissue_specificity.pdf",
       fig3_top, width = 14, height = 6, device = cairo_pdf)
ggsave("figures/fig3a_gpx4_tissue_specificity.png",
       fig3_top, width = 14, height = 6, dpi = 200)
cat("\nSaved fig3a_gpx4_tissue_specificity\n")

ggsave("figures/fig3b_trx_prot_concordance.pdf",
       fig3b, width = 10, height = 5, device = cairo_pdf)
ggsave("figures/fig3b_trx_prot_concordance.png",
       fig3b, width = 10, height = 5, dpi = 200)
cat("Saved fig3b_trx_prot_concordance\n")

ggsave("figures/fig3c_ets_vs_buffering_correlation.pdf",
       fig3c, width = 12, height = 6, device = cairo_pdf)
ggsave("figures/fig3c_ets_vs_buffering_correlation.png",
       fig3c, width = 12, height = 6, dpi = 200)
cat("Saved fig3c_ets_vs_buffering_correlation\n")

# Save stats tables
write_csv(gpx4_dat, "data/gpx4_tissue_specificity.csv")
write_csv(concordance_full, "data/trx_prot_concordance.csv")
write_csv(fig3c_dat, "data/ets_buffering_correlation_data.csv")
write_csv(reg_stats, "data/ets_buffering_regression_stats.csv")

cat("\nAll files saved.\n")
cat("\n=== KEY RESULTS SUMMARY ===\n")
cat("\n1. GPx4 tissue specificity:\n")
cat("   BAT male: logFC =", round(filter(gpx4_dat, tissue=="BAT", sex=="male")$gpx4_lfc, 3),
    " adj_p =", signif(filter(gpx4_dat, tissue=="BAT", sex=="male")$gpx4_adj_p, 3), "\n")
cat("   BLOOD male: logFC =", round(filter(gpx4_dat, tissue=="BLOOD", sex=="male")$gpx4_lfc, 3),
    " adj_p =", signif(filter(gpx4_dat, tissue=="BLOOD", sex=="male")$gpx4_adj_p, 3), "\n")
cat("   No other tissue significant at 8w\n")
cat("\n2. Transcript-protein concordance (ETS, 8w):\n")
print(concordance_full %>% select(tissue, sex_label, r_ets, p_ets) %>% arrange(sex_label, desc(r_ets)))
cat("\n3. ETS vs. buffering regression:\n")
print(reg_stats %>% select(sex_label, slope, r2, p_slope_lt1))
