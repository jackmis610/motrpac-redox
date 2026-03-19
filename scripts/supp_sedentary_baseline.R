# ============================================================
# Fig 4 — Sedentary baseline ETS vs buffering scaling
# Left panel:  baseline expression (NORM_DATA, control animals)
# Right panel: training-induced change (logFC, 8w, from fig3c)
# ============================================================

suppressPackageStartupMessages({
  library(MotrpacRatTraining6moData)
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
})

source("scripts/neufer_geneset.R")

# ── Gene category map ──────────────────────────────────────
gene_category <- imap_dfr(neufer_geneset, ~ tibble(
  gene_symbol       = .x,
  framework_category = .y
))

ets_cats <- c("ets_CI", "ets_CIII", "ets_CIV", "cyt_c", "atp_synthase")
buf_cats <- c("gsh", "trx", "nnt", "sod")

ets_genes <- gene_category %>% filter(framework_category %in% ets_cats) %>% pull(gene_symbol) %>% unique()
buf_genes <- gene_category %>% filter(framework_category %in% buf_cats) %>% pull(gene_symbol) %>% unique()

cat("ETS genes:", length(ets_genes), "| Buffering genes:", length(buf_genes), "\n")

# ── PHENO: control viallabels with sex ─────────────────────
data("PHENO")
control_pheno <- PHENO %>%
  filter(group == "control") %>%
  select(viallabel, sex, tissue) %>%
  mutate(viallabel = as.character(viallabel))

cat("Control animals:", nrow(control_pheno), "\n")

# ── FEATURE_TO_GENE: filter to ETS + buffering genes ───────
data("FEATURE_TO_GENE")
neufer_annot <- FEATURE_TO_GENE %>%
  filter(!is.na(gene_symbol),
         gene_symbol %in% c(ets_genes, buf_genes)) %>%
  select(feature_ID, gene_symbol) %>%
  distinct() %>%
  mutate(gene_class = if_else(gene_symbol %in% ets_genes, "ets", "buffering"))

cat("Annotation rows:", nrow(neufer_annot), "\n")

# ── Tissue → NORM_DATA object name map ─────────────────────
tissue_norm_map <- tribble(
  ~tissue,   ~obj_name,
  "SKM-GN",  "TRNSCRPT_SKMGN_NORM_DATA",
  "SKM-VL",  "TRNSCRPT_SKMVL_NORM_DATA",
  "HEART",   "TRNSCRPT_HEART_NORM_DATA",
  "BAT",     "TRNSCRPT_BAT_NORM_DATA",
  "WAT-SC",  "TRNSCRPT_WATSC_NORM_DATA",
  "LIVER",   "TRNSCRPT_LIVER_NORM_DATA",
  "KIDNEY",  "TRNSCRPT_KIDNEY_NORM_DATA",
  "LUNG",    "TRNSCRPT_LUNG_NORM_DATA",
  "ADRNL",   "TRNSCRPT_ADRNL_NORM_DATA",
  "BLOOD",   "TRNSCRPT_BLOOD_NORM_DATA",
  "CORTEX",  "TRNSCRPT_CORTEX_NORM_DATA",
  "HIPPOC",  "TRNSCRPT_HIPPOC_NORM_DATA",
  "HYPOTH",  "TRNSCRPT_HYPOTH_NORM_DATA",
  "VENACV",  "TRNSCRPT_VENACV_NORM_DATA",
  "COLON",   "TRNSCRPT_COLON_NORM_DATA",
  "SMLINT",  "TRNSCRPT_SMLINT_NORM_DATA",
  "SPLEEN",  "TRNSCRPT_SPLEEN_NORM_DATA",
  "OVARY",   "TRNSCRPT_OVARY_NORM_DATA",
  "TESTES",  "TRNSCRPT_TESTES_NORM_DATA"
)

# ── Extract baseline expression for each tissue ─────────────
extract_baseline <- function(tissue_label, obj_name) {
  e <- new.env()
  data(list = obj_name, package = "MotrpacRatTraining6moData", envir = e)
  norm_df <- get(obj_name, envir = e)

  # Identify sample columns (viallabels — numeric-named columns)
  meta_cols <- c("feature", "feature_ID", "tissue", "assay")
  sample_cols <- setdiff(colnames(norm_df), meta_cols)

  # Filter to control viallabels for this tissue
  ctrl_vials <- control_pheno %>%
    filter(tissue == tissue_label) %>%
    pull(viallabel)

  sample_cols_ctrl <- intersect(sample_cols, ctrl_vials)

  if (length(sample_cols_ctrl) == 0) {
    cat("  No control vials matched for", tissue_label, "— skipping\n")
    return(NULL)
  }

  # Inner join on feature_ID to get only ETS+buffering genes
  norm_filt <- norm_df %>%
    select(feature_ID, all_of(sample_cols_ctrl)) %>%
    inner_join(neufer_annot, by = "feature_ID")

  if (nrow(norm_filt) == 0) return(NULL)

  # Pivot long; join sex from control_pheno
  long <- norm_filt %>%
    pivot_longer(
      cols      = all_of(sample_cols_ctrl),
      names_to  = "viallabel",
      values_to = "expr"
    ) %>%
    left_join(control_pheno %>% filter(tissue == tissue_label) %>% select(viallabel, sex),
              by = "viallabel") %>%
    filter(!is.na(expr), !is.na(sex))

  # Mean per gene × sex (aggregate across features for same gene)
  per_gene <- long %>%
    group_by(gene_symbol, gene_class, sex) %>%
    summarise(gene_mean = mean(expr, na.rm = TRUE), .groups = "drop")

  # Mean per gene class × sex
  per_class <- per_gene %>%
    group_by(gene_class, sex) %>%
    summarise(
      class_mean = mean(gene_mean, na.rm = TRUE),
      class_se   = sd(gene_mean, na.rm = TRUE) / sqrt(n()),
      n_genes    = n(),
      .groups    = "drop"
    ) %>%
    mutate(tissue = tissue_label)

  per_class
}

cat("\nExtracting baseline expression per tissue...\n")
baseline_list <- map2(tissue_norm_map$tissue, tissue_norm_map$obj_name,
                      ~ { cat(" ", .x, "\n"); extract_baseline(.x, .y) })

baseline_long <- bind_rows(baseline_list)
cat("\nBaseline rows:", nrow(baseline_long), "\n")

# ── Pivot to wide (ets_mean vs buf_mean) ───────────────────
baseline_wide <- baseline_long %>%
  select(tissue, sex, gene_class, class_mean, class_se) %>%
  pivot_wider(
    names_from  = gene_class,
    values_from = c(class_mean, class_se)
  ) %>%
  rename(
    ets_mean = class_mean_ets,
    buf_mean = class_mean_buffering,
    ets_se   = class_se_ets,
    buf_se   = class_se_buffering
  ) %>%
  mutate(sex_label = str_to_title(sex)) %>%
  filter(!is.na(ets_mean), !is.na(buf_mean))

write_csv(baseline_wide, "data/sedentary_baseline_ets_buffering.csv")
cat("Saved data/sedentary_baseline_ets_buffering.csv\n")
print(baseline_wide %>% select(tissue, sex_label, ets_mean, buf_mean))

# ── Regression: baseline (expect slope ≈ 1) ────────────────
reg_results <- baseline_wide %>%
  group_by(sex_label) %>%
  group_modify(~ {
    fit  <- lm(buf_mean ~ ets_mean, data = .x)
    cf   <- coef(fit)
    se   <- summary(fit)$coefficients["ets_mean", "Std. Error"]
    r2   <- summary(fit)$r.squared
    n    <- nrow(.x)
    # One-tailed test: H0: slope = 1, H1: slope < 1
    t_s1 <- (cf["ets_mean"] - 1) / se
    p_s1 <- pt(t_s1, df = n - 2)
    tibble(
      intercept = cf["(Intercept)"],
      slope     = cf["ets_mean"],
      slope_se  = se,
      r_squared = r2,
      t_vs1     = t_s1,
      p_vs1     = p_s1,
      n         = n
    )
  })

cat("\n=== Baseline regression (slope vs 1) ===\n")
print(reg_results)

# ── Training-induced change data (from fig3c) ───────────────
trained_wide <- read_csv("data/ets_buffering_correlation_data.csv",
                         show_col_types = FALSE) %>%
  filter(!is.na(ets_mean), !is.na(buf_mean))

reg_trained <- trained_wide %>%
  group_by(sex_label) %>%
  group_modify(~ {
    fit <- lm(buf_mean ~ ets_mean, data = .x)
    cf  <- coef(fit)
    se  <- summary(fit)$coefficients["ets_mean", "Std. Error"]
    r2  <- summary(fit)$r.squared
    n   <- nrow(.x)
    t_s1 <- (cf["ets_mean"] - 1) / se
    p_s1 <- pt(t_s1, df = n - 2)
    tibble(
      intercept = cf["(Intercept)"],
      slope     = cf["ets_mean"],
      slope_se  = se,
      r_squared = r2,
      t_vs1     = t_s1,
      p_vs1     = p_s1,
      n         = n
    )
  })

cat("\n=== Trained regression (slope vs 1) ===\n")
print(reg_trained)

# ── Shared aesthetics ───────────────────────────────────────
SEX_COL   <- c(Male = "#2166AC", Female = "#D6604D")
SEX_SHAPE <- c(Male = 16, Female = 17)

# Helper: annotation label with slope + p-value
slope_label <- function(reg_df, sex) {
  r <- reg_df %>% filter(sex_label == sex)
  p_str <- if (r$p_vs1 < 0.001) "p<0.001" else sprintf("p=%.3f", r$p_vs1)
  sprintf("\u03b2=%.2f (%s vs 1)", r$slope, p_str)
}

# ── Panel A: Sedentary baseline ────────────────────────────
build_baseline_annotations <- function(df, reg_df) {
  df %>%
    group_by(sex_label) %>%
    group_modify(~ {
      r <- reg_df %>% filter(sex_label == .y$sex_label)
      x_rng <- range(.x$ets_mean)
      x_seq <- seq(x_rng[1], x_rng[2], length.out = 50)
      tibble(x = x_seq, y = r$intercept + r$slope * x_seq)
    })
}

anno_base <- build_baseline_annotations(baseline_wide, reg_results)

male_lbl_base   <- slope_label(reg_results, "Male")
female_lbl_base <- slope_label(reg_results, "Female")

pA <- ggplot(baseline_wide, aes(x = ets_mean, y = buf_mean,
                                 colour = sex_label, shape = sex_label)) +
  # 1:1 reference line through origin
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey60", linewidth = 0.6) +
  # Regression lines
  geom_line(data = anno_base,
            aes(x = x, y = y, colour = sex_label),
            linewidth = 0.9, show.legend = FALSE) +
  # Error bars
  geom_errorbar(aes(ymin = buf_mean - buf_se, ymax = buf_mean + buf_se),
                width = 0, alpha = 0.4, linewidth = 0.5) +
  geom_errorbar(aes(xmin = ets_mean - ets_se, xmax = ets_mean + ets_se),
                orientation = "y", width = 0, alpha = 0.4, linewidth = 0.5) +
  # Points
  geom_point(size = 2.8, alpha = 0.9) +
  # Labels
  geom_text_repel(aes(label = tissue), size = 2.4, max.overlaps = 20,
                  show.legend = FALSE) +
  # Slope annotations
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5,
           label = male_lbl_base, colour = SEX_COL["Male"],
           size = 3, fontface = "italic") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 3.2,
           label = female_lbl_base, colour = SEX_COL["Female"],
           size = 3, fontface = "italic") +
  scale_colour_manual(values = SEX_COL, name = NULL) +
  scale_shape_manual(values = SEX_SHAPE, name = NULL) +
  labs(
    title    = "Sedentary baseline",
    subtitle = "Mean normalized expression (control animals)",
    x        = "ETS index (mean norm. expr.)",
    y        = "Buffering index (mean norm. expr.)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(colour = "grey40", size = 9),
    legend.position = c(0.85, 0.15)
  )

# ── Panel B: Training-induced change (recreate fig3c style) ─
build_trained_annotations <- function(df, reg_df) {
  df %>%
    group_by(sex_label) %>%
    group_modify(~ {
      r <- reg_df %>% filter(sex_label == .y$sex_label)
      x_rng <- range(.x$ets_mean)
      x_seq <- seq(x_rng[1], x_rng[2], length.out = 50)
      tibble(x = x_seq, y = r$intercept + r$slope * x_seq)
    })
}

anno_train <- build_trained_annotations(trained_wide, reg_trained)

male_lbl_train   <- slope_label(reg_trained, "Male")
female_lbl_train <- slope_label(reg_trained, "Female")

pB <- ggplot(trained_wide, aes(x = ets_mean, y = buf_mean,
                                colour = sex_label, shape = sex_label)) +
  # 1:1 reference
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey60", linewidth = 0.6) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.4) +
  # Regression lines
  geom_line(data = anno_train,
            aes(x = x, y = y, colour = sex_label),
            linewidth = 0.9, show.legend = FALSE) +
  # Error bars
  geom_errorbar(aes(ymin = buf_mean - buf_se, ymax = buf_mean + buf_se),
                width = 0, alpha = 0.4, linewidth = 0.5) +
  geom_errorbar(aes(xmin = ets_mean - ets_se, xmax = ets_mean + ets_se),
                orientation = "y", width = 0, alpha = 0.4, linewidth = 0.5) +
  # Points
  geom_point(size = 2.8, alpha = 0.9) +
  # Labels
  geom_text_repel(aes(label = tissue), size = 2.4, max.overlaps = 20,
                  show.legend = FALSE) +
  # Slope annotations
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5,
           label = male_lbl_train, colour = SEX_COL["Male"],
           size = 3, fontface = "italic") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 3.2,
           label = female_lbl_train, colour = SEX_COL["Female"],
           size = 3, fontface = "italic") +
  scale_colour_manual(values = SEX_COL, name = NULL) +
  scale_shape_manual(values = SEX_SHAPE, name = NULL) +
  labs(
    title    = "8-week training response",
    subtitle = "Mean log\u2082FC (trained vs. sedentary, TRNSCRPT)",
    x        = "ETS index (mean log\u2082FC)",
    y        = "Buffering index (mean log\u2082FC)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(colour = "grey40", size = 9),
    legend.position = c(0.85, 0.15)
  )

# ── Assemble two-panel figure ───────────────────────────────
fig <- pA + pB +
  plot_layout(ncol = 2) +
  plot_annotation(
    title   = "Endurance training couples redox buffering to oxidative capacity across tissues, but buffering scales sub-linearly",
    caption = "Dashed line = 1:1 reference. Regression slope (\u03b2) tested vs. H\u2080: \u03b2 = 1 (one-tailed t-test).",
    theme   = theme(
      plot.title   = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.caption = element_text(colour = "grey50", size = 8, hjust = 0)
    )
  )

# ── Save ────────────────────────────────────────────────────
ggsave("figures/fig4_sedentary_vs_trained_scaling.pdf",
       fig, width = 14, height = 6.5, device = cairo_pdf)
ggsave("figures/fig4_sedentary_vs_trained_scaling.png",
       fig, width = 14, height = 6.5, dpi = 200)

cat("\nSaved figures/fig4_sedentary_vs_trained_scaling (.pdf/.png)\n")
cat("Done.\n")
