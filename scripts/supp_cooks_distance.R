# ============================================================
# Cook's distance analysis on main ETS→buffering regression
# Males, 8w, TRNSCRPT, cross-tissue
# Addresses influential-points reviewer critique
# Output: figures/fig_cooks_distance.pdf/.png
#         data/regression_diagnostics.csv
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
})

ets_cats <- c("ets_CI","ets_CIII","ets_CIV","cyt_c","atp_synthase")
buf_cats <- c("gsh","trx","nnt","sod")

# ── Load and filter ──────────────────────────────────────────
cat("Loading data...\n")
da <- read_csv("data/neufer_da_all.csv", show_col_types = FALSE)

make_indices <- function(sex_filter) {
  da %>%
    filter(assay == "TRNSCRPT",
           comparison_group == "8w",
           sex == sex_filter,
           framework_category %in% c(ets_cats, buf_cats)) %>%
    mutate(gene_class = if_else(framework_category %in% ets_cats,
                                "ets", "buffering")) %>%
    group_by(tissue, gene_symbol, gene_class) %>%
    summarise(gene_logFC = mean(logFC, na.rm = TRUE), .groups = "drop") %>%
    group_by(tissue, gene_class) %>%
    summarise(index_mean = mean(gene_logFC, na.rm = TRUE),
              index_se   = sd(gene_logFC, na.rm = TRUE) / sqrt(n()),
              n_genes    = n(), .groups = "drop") %>%
    pivot_wider(names_from  = gene_class,
                values_from = c(index_mean, index_se, n_genes)) %>%
    rename(ets_mean = index_mean_ets,
           buf_mean = index_mean_buffering,
           ets_se   = index_se_ets,
           buf_se   = index_se_buffering) %>%
    filter(!is.na(ets_mean), !is.na(buf_mean))
}

idx_m <- make_indices("male")
idx_f <- make_indices("female")

cat("Male tissues:", nrow(idx_m), "\n")
cat("Female tissues:", nrow(idx_f), "\n")

# ── Fit regressions ──────────────────────────────────────────
fit_m <- lm(buf_mean ~ ets_mean, data = idx_m)
fit_f <- lm(buf_mean ~ ets_mean, data = idx_f)

slope_m  <- coef(fit_m)[2]; r2_m <- summary(fit_m)$r.squared
slope_f  <- coef(fit_f)[2]; r2_f <- summary(fit_f)$r.squared
se_m     <- summary(fit_m)$coefficients["ets_mean","Std. Error"]
se_f     <- summary(fit_f)$coefficients["ets_mean","Std. Error"]
p_m      <- pt((slope_m - 1) / se_m, df = nrow(idx_m) - 2)
p_f      <- pt((slope_f - 1) / se_f, df = nrow(idx_f) - 2)

cat(sprintf("\nMale  8w: slope=%.4f  R²=%.4f  p(β<1)=%.2e\n", slope_m, r2_m, p_m))
cat(sprintf("Female 8w: slope=%.4f  R²=%.4f  p(β<1)=%.2e\n", slope_f, r2_f, p_f))

# ── Diagnostics ──────────────────────────────────────────────
diag_m <- idx_m %>%
  mutate(
    sex         = "Male",
    fitted      = fitted(fit_m),
    residual    = residuals(fit_m),
    std_resid   = rstandard(fit_m),
    cooks_d     = cooks.distance(fit_m),
    leverage    = hatvalues(fit_m),
    n           = nrow(idx_m),
    cook_thresh = 4 / n,
    influential = cooks_d > cook_thresh
  )

diag_f <- idx_f %>%
  mutate(
    sex         = "Female",
    fitted      = fitted(fit_f),
    residual    = residuals(fit_f),
    std_resid   = rstandard(fit_f),
    cooks_d     = cooks.distance(fit_f),
    leverage    = hatvalues(fit_f),
    n           = nrow(idx_f),
    cook_thresh = 4 / n,
    influential = cooks_d > cook_thresh
  )

diag_all <- bind_rows(diag_m, diag_f)

cat("\n=== INFLUENTIAL TISSUES (Cook's D > 4/n) ===\n")
diag_all %>%
  filter(influential) %>%
  select(sex, tissue, cooks_d, cook_thresh, ets_mean, buf_mean,
         leverage, std_resid) %>%
  mutate(across(where(is.numeric), ~round(.x, 4))) %>%
  arrange(sex, desc(cooks_d)) %>%
  print(n = 30)

write_csv(diag_all %>%
            select(sex, tissue, ets_mean, buf_mean, ets_se, buf_se,
                   fitted, residual, std_resid, cooks_d, leverage,
                   cook_thresh, influential),
          "data/regression_diagnostics.csv")
cat("\nSaved data/regression_diagnostics.csv\n")

# ── Leave-one-out slopes ────────────────────────────────────
loo_slope <- function(data, sex_label) {
  map_dfr(seq_len(nrow(data)), function(i) {
    d <- data[-i, ]
    fit <- lm(buf_mean ~ ets_mean, data = d)
    s   <- summary(fit)
    se  <- s$coefficients["ets_mean","Std. Error"]
    n   <- nrow(d)
    tibble(
      sex         = sex_label,
      tissue_left_out = data$tissue[i],
      slope       = coef(fit)[2],
      r_squared   = s$r.squared,
      p_vs1       = pt((coef(fit)[2] - 1) / se, df = n - 2)
    )
  })
}

cat("\nRunning leave-one-out...\n")
loo_m <- loo_slope(idx_m, "Male")
loo_f <- loo_slope(idx_f, "Female")
loo   <- bind_rows(loo_m, loo_f)

write_csv(loo, "data/regression_loo.csv")
cat("Saved data/regression_loo.csv\n")

cat("\n=== LOO SLOPE RANGE ===\n")
loo %>%
  group_by(sex) %>%
  summarise(
    observed_slope = if_else(sex[1]=="Male", slope_m, slope_f),
    min_slope = min(slope), max_slope = max(slope),
    min_tissue = tissue_left_out[which.min(slope)],
    max_tissue = tissue_left_out[which.max(slope)],
    all_p_remain_sig = all(p_vs1 < 0.05),
    .groups = "drop"
  ) %>%
  print()

# ── Plots ────────────────────────────────────────────────────
SEX_COL <- c(Male = "#2166AC", Female = "#D6604D")

# Panel A: Cook's D bar chart
cook_bar <- function(diag, sex_col, obs_slope, obs_r2, obs_p) {
  thresh <- 4 / nrow(diag)
  ggplot(diag %>% arrange(cooks_d),
         aes(x = reorder(tissue, cooks_d), y = cooks_d,
             fill = influential)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = thresh, linetype = "dashed",
               colour = "firebrick", linewidth = 0.7) +
    annotate("text", x = 0.6, y = thresh * 1.05,
             label = sprintf("4/n = %.3f", thresh),
             hjust = 0, vjust = 0, colour = "firebrick", size = 3) +
    scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick"),
                      guide = "none") +
    coord_flip() +
    labs(
      x     = NULL,
      y     = "Cook's distance",
      title = sprintf("%s — Cook's distance per tissue\n\u03b2=%.3f, R\u00b2=%.3f, p(\u03b2<1)=%.2e",
                      unique(diag$sex), obs_slope, obs_r2, obs_p)
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 9),
      axis.text.y = element_text(size = 8)
    )
}

pA <- cook_bar(diag_m, SEX_COL["Male"],   slope_m, r2_m, p_m)
pB <- cook_bar(diag_f, SEX_COL["Female"], slope_f, r2_f, p_f)

# Panel C: LOO slope stability — Males
loo_range_m <- range(loo_m$slope)
pc <- ggplot(loo_m, aes(x = reorder(tissue_left_out, slope), y = slope)) +
  geom_hline(yintercept = slope_m, linetype = "solid",
             colour = SEX_COL["Male"], linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  geom_point(aes(colour = p_vs1 < 0.05), size = 2.5) +
  geom_segment(aes(xend = tissue_left_out, y = slope_m, yend = slope,
                   colour = p_vs1 < 0.05), linewidth = 0.5) +
  scale_colour_manual(
    values = c("TRUE" = SEX_COL["Male"], "FALSE" = "firebrick"),
    labels = c("TRUE" = "Still sig. (p<0.05)", "FALSE" = "No longer sig."),
    name   = NULL
  ) +
  coord_flip() +
  annotate("text", x = 0.5, y = slope_m + 0.02,
           label = sprintf("Full model\n\u03b2=%.3f", slope_m),
           hjust = 0, vjust = 0, size = 2.8,
           colour = SEX_COL["Male"], fontface = "italic") +
  labs(
    title = "Males: leave-one-tissue-out slope stability",
    x = "Tissue left out", y = "Slope (\u03b2)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 9),
    legend.position = c(0.75, 0.15),
    legend.text     = element_text(size = 8),
    axis.text.y     = element_text(size = 8)
  )

# Panel D: LOO slope stability — Females
pd <- ggplot(loo_f, aes(x = reorder(tissue_left_out, slope), y = slope)) +
  geom_hline(yintercept = slope_f, linetype = "solid",
             colour = SEX_COL["Female"], linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  geom_point(aes(colour = p_vs1 < 0.05), size = 2.5) +
  geom_segment(aes(xend = tissue_left_out, y = slope_f, yend = slope,
                   colour = p_vs1 < 0.05), linewidth = 0.5) +
  scale_colour_manual(
    values = c("TRUE" = SEX_COL["Female"], "FALSE" = "firebrick"),
    labels = c("TRUE" = "Still sig.", "FALSE" = "No longer sig."),
    name   = NULL
  ) +
  coord_flip() +
  annotate("text", x = 0.5, y = slope_f + 0.015,
           label = sprintf("Full model\n\u03b2=%.3f", slope_f),
           hjust = 0, vjust = 0, size = 2.8,
           colour = SEX_COL["Female"], fontface = "italic") +
  labs(
    title = "Females: leave-one-tissue-out slope stability",
    x = NULL, y = "Slope (\u03b2)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 9),
    legend.position = c(0.75, 0.15),
    legend.text     = element_text(size = 8),
    axis.text.y     = element_text(size = 8)
  )

fig <- (pA + pB) / (pc + pd) +
  plot_annotation(
    title    = "Regression diagnostics: Cook's distance & leave-one-out stability",
    subtitle = "TRNSCRPT, 8w, cross-tissue (n=18 tissues). Dashed line = \u03b2=1 null.",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(colour = "grey40", size = 9, hjust = 0.5)
    )
  )

ggsave("figures/fig_cooks_distance.pdf",
       fig, width = 14, height = 12, device = cairo_pdf)
ggsave("figures/fig_cooks_distance.png",
       fig, width = 14, height = 12, dpi = 180)

cat("\nSaved figures/fig_cooks_distance.pdf\n")
cat("Saved figures/fig_cooks_distance.png\n")
cat("Done.\n")
