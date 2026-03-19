# ============================================================
# Fig temporal — ETS vs buffering regression at each time point
# fig_temporal_slope_males.pdf   — males only, 4-panel scatter + trajectory
# fig_temporal_slope_both.pdf    — both sexes supplementary
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
})

# ── Gene categories ─────────────────────────────────────────
ets_cats <- c("ets_CI", "ets_CIII", "ets_CIV", "cyt_c", "atp_synthase")
buf_cats <- c("gsh", "trx", "nnt", "sod")

TIMEPOINTS <- c("1w", "2w", "4w", "8w")
TP_LABELS  <- c("1w" = "1 week", "2w" = "2 weeks",
                "4w" = "4 weeks", "8w" = "8 weeks")
SEX_COL    <- c(Male = "#2166AC", Female = "#D6604D")
SEX_SHAPE  <- c(Male = 16, Female = 17)

TISSUE_ORDER <- c(
  "SKM-GN","SKM-VL","HEART","BAT","WAT-SC","LIVER","KIDNEY","LUNG",
  "ADRNL","BLOOD","CORTEX","HIPPOC","HYPOTH","VENACV",
  "COLON","SMLINT","SPLEEN","OVARY","TESTES"
)

# ── Load data ────────────────────────────────────────────────
cat("Loading data...\n")
da <- read_csv("data/neufer_da_all.csv", show_col_types = FALSE)

# ── Compute tissue-level indices per sex × timepoint ────────
make_indices <- function(data) {
  data %>%
    filter(assay == "TRNSCRPT") %>%
    mutate(
      gene_class = case_when(
        framework_category %in% ets_cats ~ "ets",
        framework_category %in% buf_cats ~ "buffering",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(gene_class)) %>%
    # gene-level mean across features
    group_by(tissue, sex, comparison_group, gene_symbol, gene_class) %>%
    summarise(gene_logFC = mean(logFC, na.rm = TRUE), .groups = "drop") %>%
    # index mean and SE across genes
    group_by(tissue, sex, comparison_group, gene_class) %>%
    summarise(
      index_mean = mean(gene_logFC, na.rm = TRUE),
      index_se   = sd(gene_logFC, na.rm = TRUE) / sqrt(n()),
      n_genes    = n(),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from  = gene_class,
      values_from = c(index_mean, index_se, n_genes)
    ) %>%
    rename(
      ets_mean = index_mean_ets,
      buf_mean = index_mean_buffering,
      ets_se   = index_se_ets,
      buf_se   = index_se_buffering
    ) %>%
    mutate(
      sex_label = str_to_title(sex),
      timepoint = factor(comparison_group, levels = TIMEPOINTS)
    ) %>%
    filter(!is.na(ets_mean), !is.na(buf_mean))
}

indices <- make_indices(da)
cat("Index rows:", nrow(indices), "\n")

# ── Regression at each timepoint × sex ──────────────────────
fit_tp_reg <- function(df) {
  df %>%
    group_by(sex_label, timepoint) %>%
    group_modify(~ {
      if (nrow(.x) < 4) return(tibble())
      fit  <- lm(buf_mean ~ ets_mean, data = .x)
      cf   <- coef(fit)
      s    <- summary(fit)
      se   <- s$coefficients["ets_mean", "Std. Error"]
      n    <- nrow(.x)
      t_s1 <- (cf["ets_mean"] - 1) / se
      p_s1 <- pt(t_s1, df = n - 2)
      ci_lo <- cf["ets_mean"] - qt(0.975, df = n - 2) * se
      ci_hi <- cf["ets_mean"] + qt(0.975, df = n - 2) * se
      tibble(
        intercept = cf["(Intercept)"],
        slope     = cf["ets_mean"],
        slope_se  = se,
        ci_lo     = ci_lo,
        ci_hi     = ci_hi,
        r_squared = s$r.squared,
        p_vs1     = p_s1,
        n         = n
      )
    })
}

reg <- fit_tp_reg(indices)
cat("\n=== Regression by timepoint × sex ===\n")
reg %>%
  select(sex_label, timepoint, slope, slope_se, r_squared, p_vs1, n) %>%
  as.data.frame() %>% print()

write_csv(reg, "data/ets_buffering_temporal_regression.csv")
cat("Saved data/ets_buffering_temporal_regression.csv\n")

# ── Helper: slope annotation label ─────────────────────────
slope_lbl <- function(r) {
  p_str <- if (r$p_vs1 < 0.001) "p<0.001" else sprintf("p=%.3f", r$p_vs1)
  sprintf("\u03b2=%.2f, R\u00b2=%.2f\n(%s vs 1)", r$slope, r$r_squared, p_str)
}

# ── Helper: regression line data ────────────────────────────
reg_lines <- function(df, reg_df) {
  df %>%
    group_by(sex_label, timepoint) %>%
    group_modify(~ {
      r <- reg_df %>%
        filter(sex_label == .y$sex_label, timepoint == .y$timepoint)
      if (nrow(r) == 0) return(tibble())
      xs <- seq(min(.x$ets_mean), max(.x$ets_mean), length.out = 50)
      tibble(x = xs, y = r$intercept + r$slope * xs)
    })
}

# ── Function: build 4-panel scatter for one sex ─────────────
build_scatter_panels <- function(sex_filter, indices, reg) {

  df_sex  <- indices %>% filter(sex_label == sex_filter)
  reg_sex <- reg    %>% filter(sex_label == sex_filter)
  lines   <- reg_lines(df_sex, reg_sex)
  col     <- SEX_COL[sex_filter]

  panels <- map(TIMEPOINTS, function(tp) {
    d  <- df_sex  %>% filter(timepoint == tp)
    r  <- reg_sex %>% filter(timepoint == tp)
    ln <- lines   %>% filter(timepoint == tp)

    if (nrow(d) == 0 || nrow(r) == 0) return(NULL)

    p_str <- if (r$p_vs1 < 0.001) "p<0.001" else sprintf("p=%.3f", r$p_vs1)
    ann   <- sprintf("\u03b2=%.2f\nR\u00b2=%.2f\n%s vs 1", r$slope, r$r_squared, p_str)

    ggplot(d, aes(x = ets_mean, y = buf_mean)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  colour = "grey60", linewidth = 0.5) +
      geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
      geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
      geom_line(data = ln, aes(x = x, y = y),
                colour = col, linewidth = 0.9, inherit.aes = FALSE) +
      geom_errorbar(aes(ymin = buf_mean - buf_se, ymax = buf_mean + buf_se),
                    width = 0, alpha = 0.35, colour = col, linewidth = 0.45) +
      geom_errorbar(aes(xmin = ets_mean - ets_se, xmax = ets_mean + ets_se),
                    orientation = "y", width = 0, alpha = 0.35,
                    colour = col, linewidth = 0.45) +
      geom_point(colour = col, size = 2.5, alpha = 0.9) +
      geom_text_repel(aes(label = tissue), size = 2.1,
                      colour = col, max.overlaps = 18,
                      segment.colour = "grey60", segment.size = 0.3) +
      annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.4,
               label = ann, size = 2.8, fontface = "italic",
               colour = col) +
      scale_x_continuous(breaks = scales::pretty_breaks(4)) +
      scale_y_continuous(breaks = scales::pretty_breaks(4)) +
      labs(
        title    = TP_LABELS[tp],
        x        = "ETS index (mean log\u2082FC)",
        y        = "Buffering index (mean log\u2082FC)"
      ) +
      theme_classic(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
        axis.title = element_text(size = 9)
      )
  })
  panels[!sapply(panels, is.null)]
}

# ── Function: trajectory plot (slope + R² over time) ────────
build_trajectory <- function(reg, sex_filter = NULL) {

  df <- if (!is.null(sex_filter)) {
    reg %>% filter(sex_label %in% sex_filter)
  } else {
    reg
  }

  df <- df %>%
    mutate(tp_num = as.integer(timepoint) * 2)  # 1w=2, 2w=4, 4w=8w... use factor index

  # Map timepoint factor levels to numeric weeks
  tp_weeks <- c("1w" = 1, "2w" = 2, "4w" = 4, "8w" = 8)
  df <- df %>%
    mutate(weeks = tp_weeks[as.character(timepoint)])

  # Panel A: slope with CI ribbon
  pA <- ggplot(df, aes(x = weeks, y = slope, colour = sex_label,
                        fill = sex_label, shape = sex_label)) +
    geom_hline(yintercept = 1, linetype = "dashed",
               colour = "grey50", linewidth = 0.6) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15,
                colour = NA, show.legend = FALSE) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 3.5) +
    scale_colour_manual(values = SEX_COL, name = NULL) +
    scale_fill_manual(values   = SEX_COL, name = NULL) +
    scale_shape_manual(values  = SEX_SHAPE, name = NULL) +
    scale_x_continuous(breaks = c(1, 2, 4, 8),
                       labels = c("1w", "2w", "4w", "8w")) +
    annotate("text", x = 8.2, y = 1.02, hjust = 0, vjust = 0,
             label = "1:1", colour = "grey40", size = 3) +
    labs(
      title = "Regression slope (\u03b2) over training",
      x     = NULL,
      y     = "Slope (\u03b2, ETS \u2192 buffering)"
    ) +
    coord_cartesian(xlim = c(0.7, 9)) +
    theme_classic(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 11),
      legend.position = c(0.15, 0.85),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank()
    )

  # Panel B: R²
  pB <- ggplot(df, aes(x = weeks, y = r_squared, colour = sex_label,
                        shape = sex_label)) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 3.5) +
    scale_colour_manual(values = SEX_COL, name = NULL) +
    scale_shape_manual(values  = SEX_SHAPE, name = NULL) +
    scale_x_continuous(breaks = c(1, 2, 4, 8),
                       labels = c("1w", "2w", "4w", "8w")) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(
      title = "Variance explained (R\u00b2)",
      x     = "Training duration",
      y     = "R\u00b2"
    ) +
    coord_cartesian(xlim = c(0.7, 9)) +
    theme_classic(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 11),
      legend.position = "none"
    )

  pA / pB + plot_layout(heights = c(2, 1))
}

# ════════════════════════════════════════════════════════════
# Figure 1: Males only
# ════════════════════════════════════════════════════════════
cat("\nBuilding males-only figure...\n")

scatter_m <- build_scatter_panels("Male", indices, reg)
traj_m    <- build_trajectory(reg, sex_filter = "Male")

# Combine using wrap_plots to avoid nested patchwork annotation issue
scatters_m <- wrap_plots(scatter_m, ncol = 2)

fig_males <- wrap_plots(
  list(scatters_m, traj_m),
  ncol = 1,
  heights = c(2, 1.2)
) +
  plot_annotation(
    title    = "ETS vs. redox buffering scaling across training — Males",
    subtitle = "TRNSCRPT logFC (19 tissues). Dashed line = 1:1 reference.",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(colour = "grey40", size = 9, hjust = 0.5)
    )
  )

ggsave("figures/fig_temporal_slope_males.pdf",
       fig_males, width = 12, height = 14, device = cairo_pdf)
ggsave("figures/fig_temporal_slope_males.png",
       fig_males, width = 12, height = 14, dpi = 180)
cat("Saved figures/fig_temporal_slope_males\n")

# ════════════════════════════════════════════════════════════
# Figure 2: Both sexes supplementary
# ════════════════════════════════════════════════════════════
cat("Building both-sexes supplementary figure...\n")

scatter_f <- build_scatter_panels("Female", indices, reg)
traj_both <- build_trajectory(reg)  # both sexes on one plot

# Interleave male/female panels: row 1 = 1w M+F, row 2 = 2w M+F, etc.
# Simpler: two columns stacked — wrap_plots handles it
col_m <- wrap_plots(scatter_m, ncol = 1)
col_f <- wrap_plots(scatter_f, ncol = 1)
scatter_grid <- wrap_plots(list(col_m, col_f), ncol = 2)

fig_both <- wrap_plots(
  list(scatter_grid, traj_both),
  ncol = 1,
  heights = c(3, 1.2)
) +
  plot_annotation(
    title    = "ETS vs. redox buffering scaling across training — Males and Females",
    subtitle = "TRNSCRPT logFC (19 tissues). Dashed line = 1:1 reference.",
    caption  = paste0(
      "\u03b2 at 8w: Males=0.64 (R\u00b2=0.75, p<0.001 vs \u03b2=1); ",
      "Females=0.14 (R\u00b2=0.32, p<0.001 vs \u03b2=1)."
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(colour = "grey40", size = 9,  hjust = 0.5),
      plot.caption  = element_text(colour = "grey50", size = 8,  hjust = 0)
    )
  )

ggsave("figures/fig_temporal_slope_both.pdf",
       fig_both, width = 18, height = 18, device = cairo_pdf)
ggsave("figures/fig_temporal_slope_both.png",
       fig_both, width = 18, height = 18, dpi = 160)
cat("Saved figures/fig_temporal_slope_both\n")

# ── Print regression table ───────────────────────────────────
cat("\n=== Full regression table ===\n")
reg %>%
  mutate(across(c(slope, slope_se, ci_lo, ci_hi, r_squared), ~round(.x, 3)),
         p_vs1 = signif(p_vs1, 3)) %>%
  as.data.frame() %>% print()

cat("\nDone.\n")
