# ============================================================
# MitoCarta 3.0 comparison figure
#
# Shows that a localization-based (MitoCarta) approach fails to
# isolate the ETS→buffering sub-circuit, while our curated
# thermodynamic framework does.
#
# Three approaches compared (males, 8w, TRNSCRPT):
#   A) Neufer framework:  ETS (32 genes) vs buffering (14 genes)
#   B) MitoCarta OXPHOS vs MitoCarta ROS/redox subset
#   C) MitoCarta OXPHOS vs ALL other MitoCarta (diluted)
#
# Output: figures/fig_mitocarta_comparison.pdf/.png
# ============================================================

suppressPackageStartupMessages({
  library(MotrpacRatTraining6moData)
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
})

# ── MitoCarta 3.0 gene sets (hard-coded from published list) ──
# Source: Rath et al. (2021) MitoCarta3.0, PMID 33174596
# MitoPathways hierarchy: OXPHOS → complexes I-V
# Antioxidants/ROS → SOD, GPx, TXN, PRDX, GSR, NNT, etc.

mitocarta_oxphos <- c(
  # Complex I (NADH:ubiquinone oxidoreductase)
  "Ndufs1","Ndufs2","Ndufs3","Ndufs4","Ndufs5","Ndufs6","Ndufs7","Ndufs8",
  "Ndufv1","Ndufv2","Ndufv3",
  "Ndufa1","Ndufa2","Ndufa3","Ndufa4","Ndufa5","Ndufa6","Ndufa7","Ndufa8",
  "Ndufa9","Ndufa10","Ndufa11","Ndufa12","Ndufa13",
  "Ndufb1","Ndufb2","Ndufb3","Ndufb4","Ndufb5","Ndufb6","Ndufb7","Ndufb8",
  "Ndufb9","Ndufb10","Ndufb11",
  "Ndufc1","Ndufc2",
  "Ndufab1",
  # Complex II (succinate:ubiquinone oxidoreductase)
  "Sdha","Sdhb","Sdhc","Sdhd",
  # Complex III (ubiquinol:cytochrome c oxidoreductase)
  "Uqcrc1","Uqcrc2","Uqcrfs1","Uqcrb","Uqcrh","Uqcrq","Uqcr10","Uqcr11",
  "Cyc1","Cycs",
  # Complex IV (cytochrome c oxidase)
  "Cox4i1","Cox4i2","Cox5a","Cox5b",
  "Cox6a1","Cox6a2","Cox6b1","Cox6b2","Cox6c",
  "Cox7a1","Cox7a2","Cox7b","Cox7c",
  "Cox8a","Cox8c",
  "Cox14","Cox16","Cox17","Cox18","Cox19","Cox20",
  # Complex V (ATP synthase)
  "Atp5f1a","Atp5f1b","Atp5f1c","Atp5f1d","Atp5f1e",
  "Atp5pb","Atp5mc1","Atp5mc2","Atp5mc3","Atp5po",
  "Atp5if1","Atp5md","Atp5mg","Atp5mpl","Atp5me"
)

# MitoCarta "Reactive oxygen species" pathway
mitocarta_ros <- c(
  "Sod1","Sod2","Sod3",
  "Gpx1","Gpx2","Gpx4","Gpx6","Gpx7","Gpx8",
  "Prdx1","Prdx2","Prdx3","Prdx4","Prdx5","Prdx6",
  "Txn1","Txn2","Txnrd1","Txnrd2","Txnrd3",
  "Gsr","Glrx","Glrx2","Glrx3","Glrx5",
  "Nnt","Gclc","Gclm","Gss","Gstp1","Gstp2",
  "Hmox1","Hmox2","Nqo1","Sqstm1",
  "Cat","Idh2","Me2","Me3"
)

# Remove OXPHOS genes from the "all other MitoCarta" pool (no overlap)
# We also need a broad "other mito" set — use a representative list
# of non-OXPHOS MitoCarta genes (metabolic, dynamics, translation, etc.)
mitocarta_other <- c(
  # TCA cycle
  "Cs","Aco2","Idh3a","Idh3b","Idh3g","Ogdh","Sucla2","Suclg1","Sdha",
  "Fh","Mdh2","Dlst","Dld","Ogdhl",
  # Beta-oxidation
  "Acadm","Acadl","Acadvl","Hadha","Hadhb","Cpt1a","Cpt1b","Cpt2",
  "Etfa","Etfb","Etfdh","Echs1","Hadh",
  # Import machinery (TIM/TOM)
  "Tomm20","Tomm22","Tomm40","Tomm70","Timm23","Timm44","Timm50",
  "Hspd1","Hspe1",
  # Dynamics
  "Mfn1","Mfn2","Opa1","Dnm1l","Fis1","Mff",
  # Translation
  "Mtif2","Tufm","Mrps22","Mrpl12","Mrpl37",
  # Metabolite carriers (SLC25)
  "Slc25a3","Slc25a4","Slc25a5","Slc25a6","Slc25a11","Slc25a12","Slc25a13",
  "Slc25a20","Slc25a27",
  # PDH complex
  "Pdha1","Pdhb","Dlat","Dld","Pdhx","Pdk1","Pdk2","Pdk3","Pdk4",
  # mtUPR / quality control
  "Lonp1","Clpp","Clpx","Yme1l1","Paraplegin",
  # Misc
  "Nnt","Idh2","Ppif","Vdac1","Vdac2","Vdac3","Mcu","Mcur1"
)
# Remove any accidental overlap with OXPHOS
mitocarta_other <- setdiff(mitocarta_other, mitocarta_oxphos)

cat("MitoCarta OXPHOS genes defined:", length(mitocarta_oxphos), "\n")
cat("MitoCarta ROS genes defined:",    length(mitocarta_ros), "\n")
cat("MitoCarta 'other' genes defined:", length(mitocarta_other), "\n")

# ── Load data and build full gene×tissue matrix ───────────────
cat("\nLoading transcriptome data...\n")
da_neufer <- read_csv("data/neufer_da_all.csv", show_col_types = FALSE)

# Also need full matrix for MitoCarta genes not in neufer set
# Load from the permutation test output matrix if available,
# otherwise rebuild for relevant genes only
all_mc_genes <- unique(c(mitocarta_oxphos, mitocarta_ros, mitocarta_other))
cat("Total unique MitoCarta genes needed:", length(all_mc_genes), "\n")

data(FEATURE_TO_GENE)
f2g <- FEATURE_TO_GENE %>%
  filter(!is.na(gene_symbol), gene_symbol %in% all_mc_genes) %>%
  select(feature_ID, gene_symbol) %>%
  distinct()

cat("MitoCarta genes found in FEATURE_TO_GENE:",
    n_distinct(f2g$gene_symbol), "\n")

# Load all TRNSCRPT tables for 8w male
all_items   <- data(package = "MotrpacRatTraining6moData")$results[, "Item"]
trnscrpt_da <- all_items[grepl("^TRNSCRPT_.*_DA$", all_items)]

cat("Loading", length(trnscrpt_da), "TRNSCRPT tables...\n")
mc_logfc <- map_dfr(trnscrpt_da, function(tn) {
  e <- new.env()
  tryCatch({
    data(list = tn, package = "MotrpacRatTraining6moData", envir = e)
    get(tn, envir = e) %>%
      filter(sex == "male", comparison_group == "8w") %>%
      select(feature_ID, tissue, logFC) %>%
      inner_join(f2g, by = "feature_ID") %>%
      group_by(tissue, gene_symbol) %>%
      summarise(logFC = mean(logFC, na.rm = TRUE), .groups = "drop")
  }, error = function(e) tibble())
})

cat("MitoCarta rows loaded:", nrow(mc_logfc), "\n")

# ── Compute tissue-level mean indices ─────────────────────────
compute_index <- function(logfc_data, gene_set, label) {
  genes_present <- intersect(gene_set, logfc_data$gene_symbol)
  logfc_data %>%
    filter(gene_symbol %in% genes_present) %>%
    group_by(tissue) %>%
    summarise(!!label := mean(logFC, na.rm = TRUE), .groups = "drop")
}

# Neufer indices (from existing neufer_da_all.csv)
ets_cats <- c("ets_CI","ets_CIII","ets_CIV","cyt_c","atp_synthase")
buf_cats <- c("gsh","trx","nnt","sod")

neufer_ets <- da_neufer %>%
  filter(assay == "TRNSCRPT", sex == "male",
         comparison_group == "8w",
         framework_category %in% ets_cats) %>%
  group_by(tissue, gene_symbol) %>%
  summarise(logFC = mean(logFC, na.rm=TRUE), .groups="drop") %>%
  group_by(tissue) %>%
  summarise(neufer_ets = mean(logFC, na.rm=TRUE), .groups="drop")

neufer_buf <- da_neufer %>%
  filter(assay == "TRNSCRPT", sex == "male",
         comparison_group == "8w",
         framework_category %in% buf_cats) %>%
  group_by(tissue, gene_symbol) %>%
  summarise(logFC = mean(logFC, na.rm=TRUE), .groups="drop") %>%
  group_by(tissue) %>%
  summarise(neufer_buf = mean(logFC, na.rm=TRUE), .groups="drop")

# MitoCarta indices
mc_oxphos   <- compute_index(mc_logfc, mitocarta_oxphos, "mc_ets")
mc_ros      <- compute_index(mc_logfc, mitocarta_ros,    "mc_ros")
mc_other    <- compute_index(mc_logfc, mitocarta_other,  "mc_other")

# Merged datasets
df_neufer <- neufer_ets %>%
  inner_join(neufer_buf, by = "tissue") %>%
  rename(ets_mean = neufer_ets, buf_mean = neufer_buf) %>%
  mutate(approach = "Neufer framework\n(thermodynamic circuit)")

df_mc_ros <- mc_oxphos %>%
  inner_join(mc_ros, by = "tissue") %>%
  rename(ets_mean = mc_ets, buf_mean = mc_ros) %>%
  mutate(approach = "MitoCarta\nOXPHOS vs ROS/redox")

df_mc_all <- mc_oxphos %>%
  inner_join(mc_other, by = "tissue") %>%
  rename(ets_mean = mc_ets, buf_mean = mc_other) %>%
  mutate(approach = "MitoCarta\nOXPHOS vs all other mito")

cat("\nGenes contributing:\n")
cat("  Neufer ETS:", n_distinct(
  da_neufer %>% filter(framework_category %in% ets_cats) %>% pull(gene_symbol)), "\n")
cat("  Neufer buf:", n_distinct(
  da_neufer %>% filter(framework_category %in% buf_cats) %>% pull(gene_symbol)), "\n")
cat("  MitoCarta OXPHOS in data:",
    n_distinct(mc_logfc %>% filter(gene_symbol %in% mitocarta_oxphos) %>% pull(gene_symbol)), "\n")
cat("  MitoCarta ROS in data:",
    n_distinct(mc_logfc %>% filter(gene_symbol %in% mitocarta_ros) %>% pull(gene_symbol)), "\n")
cat("  MitoCarta other in data:",
    n_distinct(mc_logfc %>% filter(gene_symbol %in% mitocarta_other) %>% pull(gene_symbol)), "\n")

# ── Fit regressions and extract stats ────────────────────────
fit_and_stats <- function(df) {
  fit <- lm(buf_mean ~ ets_mean, data = df)
  s   <- summary(fit)
  se  <- s$coefficients["ets_mean","Std. Error"]
  n   <- nrow(df)
  p1  <- pt((coef(fit)[2] - 1) / se, df = n - 2)
  list(fit = fit, slope = coef(fit)[2], r2 = s$r.squared,
       p_vs1 = p1, n = n, intercept = coef(fit)[1])
}

stats_n   <- fit_and_stats(df_neufer)
stats_ros <- fit_and_stats(df_mc_ros)
stats_all <- fit_and_stats(df_mc_all)

cat("\n=== COMPARISON TABLE ===\n")
comp_tbl <- tibble(
  Approach      = c("Neufer framework", "MitoCarta OXPHOS vs ROS",
                    "MitoCarta OXPHOS vs all other"),
  n_tissues     = c(stats_n$n,   stats_ros$n,   stats_all$n),
  slope         = c(stats_n$slope, stats_ros$slope, stats_all$slope),
  R2            = c(stats_n$r2,    stats_ros$r2,    stats_all$r2),
  p_slope_lt1   = c(stats_n$p_vs1, stats_ros$p_vs1, stats_all$p_vs1)
)
comp_tbl %>%
  mutate(across(c(slope, R2), ~round(.x, 3)),
         p_slope_lt1 = signif(p_slope_lt1, 3)) %>%
  print()

write_csv(comp_tbl, "data/mitocarta_comparison_stats.csv")
cat("Saved data/mitocarta_comparison_stats.csv\n")

# ── Build scatter panels ───────────────────────────────────────
APPROACH_COL <- c(
  "Neufer framework\n(thermodynamic circuit)" = "#2166AC",
  "MitoCarta\nOXPHOS vs ROS/redox"            = "#4DAC26",
  "MitoCarta\nOXPHOS vs all other mito"        = "#F1A340"
)

scatter_panel <- function(df, stats, col) {
  approach_label <- unique(df$approach)
  ann <- sprintf("\u03b2 = %.3f\nR\u00b2 = %.3f\np(\u03b2<1) = %s",
                 stats$slope, stats$r2,
                 if (stats$p_vs1 < 0.001) "<0.001"
                 else sprintf("%.3f", stats$p_vs1))

  xs <- seq(min(df$ets_mean, na.rm=TRUE),
            max(df$ets_mean, na.rm=TRUE), length.out = 50)
  line_df <- tibble(x = xs, y = stats$intercept + stats$slope * xs)

  ggplot(df, aes(x = ets_mean, y = buf_mean)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey60", linewidth = 0.5) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_line(data = line_df, aes(x = x, y = y),
              colour = col, linewidth = 1.0, inherit.aes = FALSE) +
    geom_point(colour = col, size = 2.5, alpha = 0.85) +
    geom_text_repel(aes(label = tissue), size = 2.2,
                    colour = col, max.overlaps = 16,
                    segment.colour = "grey65", segment.size = 0.3) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
             label = ann, size = 3.0, fontface = "italic",
             colour = col) +
    labs(
      title = approach_label,
      x     = "\"ETS\" index (mean log\u2082FC)",
      y     = "\"Buffering\" index (mean log\u2082FC)"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5)
    )
}

p1 <- scatter_panel(df_neufer,  stats_n,   APPROACH_COL[1])
p2 <- scatter_panel(df_mc_ros,  stats_ros, APPROACH_COL[2])
p3 <- scatter_panel(df_mc_all,  stats_all, APPROACH_COL[3])

# Panel D: R² comparison bar chart
r2_df <- tibble(
  Approach = c("Neufer\nframework", "MitoCarta\nOXPHOS\nvs ROS",
               "MitoCarta\nOXPHOS vs\nall other"),
  R2       = c(stats_n$r2, stats_ros$r2, stats_all$r2),
  slope    = c(stats_n$slope, stats_ros$slope, stats_all$slope),
  col      = unname(APPROACH_COL)
)

p4 <- ggplot(r2_df, aes(x = Approach, y = R2, fill = Approach)) +
  geom_col(width = 0.55, show.legend = FALSE) +
  geom_text(aes(label = sprintf("R\u00b2=%.3f\n\u03b2=%.3f", R2, slope)),
            vjust = -0.3, size = 3.0, fontface = "bold") +
  scale_fill_manual(values = setNames(r2_df$col, r2_df$Approach)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Cross-tissue ETS\u2013buffering coupling (R\u00b2)",
    subtitle = "Higher R\u00b2 = stronger coupling; all males, 8w TRNSCRPT",
    x = NULL, y = "R\u00b2"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title    = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(colour = "grey50", size = 8),
    axis.text.x   = element_text(size = 8, lineheight = 1.1)
  )

fig <- (p1 + p2 + p3) / (plot_spacer() + p4 + plot_spacer()) +
  plot_layout(heights = c(3, 2)) +
  plot_annotation(
    title    = "MitoCarta vs. Neufer framework: cross-tissue ETS\u2013buffering regression (Males, 8w)",
    subtitle = paste0("MitoCarta localization-based grouping vs. thermodynamic circuit-based grouping. ",
                      "Dashed line = 1:1 reference."),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(colour = "grey40", size = 9, hjust = 0.5)
    )
  )

ggsave("figures/fig_mitocarta_comparison.pdf",
       fig, width = 16, height = 13, device = cairo_pdf)
ggsave("figures/fig_mitocarta_comparison.png",
       fig, width = 16, height = 13, dpi = 180)

cat("\nSaved figures/fig_mitocarta_comparison.pdf\n")
cat("Saved figures/fig_mitocarta_comparison.png\n")
cat("Done.\n")
