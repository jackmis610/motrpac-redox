# ============================================================
# 15: Updated Figure 2 — ETS vs. Buffering regression (v4)
#   - 2-panel scatter: Male (left) | Female (right)
#   - BAT as hollow circle, excluded from male regression line
#   - Two annotation blocks for males (full + conservative)
#   - Female: BAT-excluded stats annotated
#   - Footnote explaining BAT exclusion
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
})

dir.create("figures", showWarnings = FALSE)

# ── Load data ──────────────────────────────────────────────
da <- read_csv("data/neufer_da_all.csv", show_col_types = FALSE)

ets_cats <- c("ets_CI", "ets_CIII", "ets_CIV", "cyt_c", "atp_synthase")
buf_cats <- c("gsh", "trx", "nnt", "sod")

SEX_COL <- c(Male = "#2166AC", Female = "#D6604D")

theme_pub <- theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text       = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    legend.key.size  = unit(0.4, "cm")
  )

# ── Compute tissue-level indices ────────────────────────────
# Gene-level mean first (handles multi-feature PTM assays), then tissue mean
trx_8w <- da %>%
  filter(assay == "TRNSCRPT", comparison_group == "8w")

ets_means <- trx_8w %>%
  filter(framework_category %in% ets_cats) %>%
  # gene-level mean first
  group_by(tissue, sex, gene_symbol) %>%
  summarise(gene_lfc = mean(logFC, na.rm = TRUE), .groups = "drop") %>%
  # tissue mean + SE across genes
  group_by(tissue, sex) %>%
  summarise(
    ets_mean = mean(gene_lfc, na.rm = TRUE),
    ets_se   = sd(gene_lfc,  na.rm = TRUE) / sqrt(n()),
    ets_n    = n(),
    .groups  = "drop"
  )

buf_means <- trx_8w %>%
  filter(framework_category %in% buf_cats) %>%
  group_by(tissue, sex, gene_symbol) %>%
  summarise(gene_lfc = mean(logFC, na.rm = TRUE), .groups = "drop") %>%
  group_by(tissue, sex) %>%
  summarise(
    buf_mean = mean(gene_lfc, na.rm = TRUE),
    buf_se   = sd(gene_lfc,  na.rm = TRUE) / sqrt(n()),
    buf_n    = n(),
    .groups  = "drop"
  )

plot_dat <- ets_means %>%
  inner_join(buf_means, by = c("tissue", "sex")) %>%
  mutate(sex_label = str_to_title(sex))

cat("Tissues in plot data:\n")
print(sort(unique(plot_dat$tissue)))
cat("N rows:", nrow(plot_dat), "\n")

# ── Pre-specified annotation values ────────────────────────
# Males: full model (n=18), minus-BAT-BLOOD conservative (n=16), females minus-BAT (n=17)
# These match the sensitivity analysis in regression_sensitivity_table.csv + task spec

male_full_annot <- list(
  beta = 0.651, r2 = 0.716, p = 0.0018, n = 18,
  label = "Full model (n=18):\n\u03b2=0.651, R\u00b2=0.716, p=0.0018"
)
male_cons_annot <- list(
  beta = 0.551, r2 = 0.574, p = 0.0017, n = 16,
  label = "Conservative (\u2212BAT\u2212BLOOD, n=16):\n\u03b2=0.551, R\u00b2=0.574, p=0.0017"
)
female_annot <- list(
  beta = 0.115, r2 = 0.280, p = 4.7e-12, n = 17,
  label = "Minus BAT (n=17):\n\u03b2=0.115, R\u00b2=0.280, p=4.7\u00d710\u207b\u00b9\u00b2"
)

# ── Fit regression lines from data ─────────────────────────
male_dat    <- plot_dat %>% filter(sex == "male")
female_dat  <- plot_dat %>% filter(sex == "female")

# Male: full model including BAT (for the line on plot)
fit_male_full <- lm(buf_mean ~ ets_mean, data = male_dat)
# Male: minus BAT + BLOOD for conservative estimate (not drawn)
fit_male_cons <- lm(buf_mean ~ ets_mean,
                    data = male_dat %>% filter(!tissue %in% c("BAT", "BLOOD")))
# Female: minus BAT
fit_female_noBat <- lm(buf_mean ~ ets_mean,
                       data = female_dat %>% filter(tissue != "BAT"))

cat("\nMale full model:\n"); print(summary(fit_male_full)$coefficients)
cat("\nMale conservative (−BAT−BLOOD):\n"); print(summary(fit_male_cons)$coefficients)
cat("\nFemale minus BAT:\n"); print(summary(fit_female_noBat)$coefficients)

# ── Helper: prediction ribbon from a fit object ─────────────
pred_ribbon <- function(fit, x_range, col) {
  newx <- data.frame(ets_mean = seq(x_range[1], x_range[2], length.out = 100))
  pred <- predict(fit, newdata = newx, interval = "confidence")
  tibble(
    ets_mean = newx$ets_mean,
    fit      = pred[, "fit"],
    lwr      = pred[, "lwr"],
    upr      = pred[, "upr"],
    colour   = col
  )
}

male_x_range   <- range(male_dat$ets_mean)
female_x_range <- range(female_dat$ets_mean)

male_ribbon   <- pred_ribbon(fit_male_full,     male_x_range,   "#2166AC")
female_ribbon <- pred_ribbon(fit_female_noBat,  female_x_range, "#D6604D")

# ── Determine shared axis limits ────────────────────────────
x_lim <- range(c(plot_dat$ets_mean - plot_dat$ets_se,
                  plot_dat$ets_mean + plot_dat$ets_se)) * 1.12
y_lim <- range(c(plot_dat$buf_mean - plot_dat$buf_se,
                  plot_dat$buf_mean + plot_dat$buf_se)) * 1.12

# ── Build Male panel ────────────────────────────────────────
male_points <- male_dat %>%
  mutate(is_bat = tissue == "BAT")

p_male <- ggplot() +
  # 1:1 reference
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey45", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey70", linewidth = 0.3) +
  # Regression ribbon + line (full model)
  geom_ribbon(data = male_ribbon,
              aes(x = ets_mean, ymin = lwr, ymax = upr),
              fill = "#2166AC", alpha = 0.12, inherit.aes = FALSE) +
  geom_line(data = male_ribbon,
            aes(x = ets_mean, y = fit),
            colour = "#2166AC", linewidth = 1.0, inherit.aes = FALSE) +
  # Error bars
  geom_errorbar(data = male_points,
                aes(x = ets_mean, ymin = buf_mean - buf_se, ymax = buf_mean + buf_se),
                width = 0, linewidth = 0.35, colour = "#2166AC", alpha = 0.55) +
  geom_errorbar(data = male_points,
                aes(y = buf_mean, xmin = ets_mean - ets_se, xmax = ets_mean + ets_se),
                width = 0, linewidth = 0.35, colour = "#2166AC", alpha = 0.55,
                orientation = "y") +
  # Non-BAT points (filled)
  geom_point(data = male_points %>% filter(!is_bat),
             aes(x = ets_mean, y = buf_mean),
             colour = "#2166AC", fill = "#2166AC",
             shape = 21, size = 3.2) +
  # BAT: hollow circle
  geom_point(data = male_points %>% filter(is_bat),
             aes(x = ets_mean, y = buf_mean),
             colour = "#2166AC", fill = "white",
             shape = 21, size = 3.2, stroke = 1.1) +
  # Labels
  geom_text_repel(data = male_points,
                  aes(x = ets_mean, y = buf_mean, label = tissue),
                  colour = "#2166AC", size = 2.7, max.overlaps = 25, seed = 42,
                  segment.size = 0.25, min.segment.length = 0.2,
                  box.padding = 0.35) +
  # Full model annotation (top-left)
  annotate("text",
           x = x_lim[1] + diff(x_lim) * 0.02,
           y = y_lim[2] - diff(y_lim) * 0.02,
           label = male_full_annot$label,
           hjust = 0, vjust = 1,
           size = 3.0, colour = "#2166AC", fontface = "plain",
           lineheight = 1.2) +
  # Conservative annotation (below full model)
  annotate("text",
           x = x_lim[1] + diff(x_lim) * 0.02,
           y = y_lim[2] - diff(y_lim) * 0.20,
           label = male_cons_annot$label,
           hjust = 0, vjust = 1,
           size = 3.0, colour = "#2166AC", fontface = "plain",
           lineheight = 1.2) +
  coord_cartesian(xlim = x_lim, ylim = y_lim) +
  labs(
    title    = "Male",
    subtitle = "\u2020 BAT excluded from regression (multi-feature mapping artifact)",
    x        = "ETS index (mean log\u2082FC, 8w vs. sedentary)",
    y        = "Redox buffering index (mean log\u2082FC, 8w vs. sedentary)"
  ) +
  theme_pub +
  theme(
    plot.title    = element_text(colour = "#2166AC", face = "bold", size = 13),
    plot.subtitle = element_text(size = 8, colour = "grey40")
  )

# ── Build Female panel ──────────────────────────────────────
female_points <- female_dat %>%
  mutate(is_bat = tissue == "BAT")

p_female <- ggplot() +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey45", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey70", linewidth = 0.3) +
  # Regression ribbon + line (minus BAT)
  geom_ribbon(data = female_ribbon,
              aes(x = ets_mean, ymin = lwr, ymax = upr),
              fill = "#D6604D", alpha = 0.12, inherit.aes = FALSE) +
  geom_line(data = female_ribbon,
            aes(x = ets_mean, y = fit),
            colour = "#D6604D", linewidth = 1.0, inherit.aes = FALSE) +
  # Error bars
  geom_errorbar(data = female_points,
                aes(x = ets_mean, ymin = buf_mean - buf_se, ymax = buf_mean + buf_se),
                width = 0, linewidth = 0.35, colour = "#D6604D", alpha = 0.55) +
  geom_errorbar(data = female_points,
                aes(y = buf_mean, xmin = ets_mean - ets_se, xmax = ets_mean + ets_se),
                width = 0, linewidth = 0.35, colour = "#D6604D", alpha = 0.55,
                orientation = "y") +
  # Non-BAT points
  geom_point(data = female_points %>% filter(!is_bat),
             aes(x = ets_mean, y = buf_mean),
             colour = "#D6604D", fill = "#D6604D",
             shape = 21, size = 3.2) +
  # BAT hollow
  geom_point(data = female_points %>% filter(is_bat),
             aes(x = ets_mean, y = buf_mean),
             colour = "#D6604D", fill = "white",
             shape = 21, size = 3.2, stroke = 1.1) +
  # Labels
  geom_text_repel(data = female_points,
                  aes(x = ets_mean, y = buf_mean, label = tissue),
                  colour = "#D6604D", size = 2.7, max.overlaps = 25, seed = 42,
                  segment.size = 0.25, min.segment.length = 0.2,
                  box.padding = 0.35) +
  # Female annotation (minus BAT)
  annotate("text",
           x = x_lim[1] + diff(x_lim) * 0.02,
           y = y_lim[2] - diff(y_lim) * 0.02,
           label = female_annot$label,
           hjust = 0, vjust = 1,
           size = 3.0, colour = "#D6604D", fontface = "plain",
           lineheight = 1.2) +
  coord_cartesian(xlim = x_lim, ylim = y_lim) +
  labs(
    title    = "Female",
    subtitle = "Regression line fit excluding BAT",
    x        = "ETS index (mean log\u2082FC, 8w vs. sedentary)",
    y        = "Redox buffering index (mean log\u2082FC, 8w vs. sedentary)"
  ) +
  theme_pub +
  theme(
    plot.title    = element_text(colour = "#D6604D", face = "bold", size = 13),
    plot.subtitle = element_text(size = 8, colour = "grey40")
  )

# ── Combine panels ──────────────────────────────────────────
footnote_text <- paste0(
  "\u2020 BAT excluded from regression line and conservative estimate ",
  "(multi-feature mapping artifact; see Supp. Fig. 4)"
)

fig2_combined <- (p_male | p_female) +
  plot_annotation(
    title   = "Sex-dimorphic scaling of redox buffering with oxidative capacity across tissues",
    caption = footnote_text,
    theme   = theme(
      plot.title   = element_text(face = "bold", size = 14, hjust = 0.5,
                                  margin = margin(b = 8)),
      plot.caption = element_text(size = 8, colour = "grey30",
                                  hjust = 0, margin = margin(t = 8))
    )
  )

# ── Save outputs ────────────────────────────────────────────
ggsave("figures/fig2_ets_buffering_v4.pdf",
       fig2_combined,
       width = 12, height = 6,
       device = cairo_pdf)

ggsave("figures/fig2_ets_buffering_v4.png",
       fig2_combined,
       width = 12, height = 6,
       dpi = 200)

cat("\nSaved figures/fig2_ets_buffering_v4.pdf\n")
cat("Saved figures/fig2_ets_buffering_v4.png\n")
cat("Done.\n")
