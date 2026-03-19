# =============================================================================
# 07_human_data.R
# Purpose : Cross-reference rat findings with human MoTrPAC acute exercise data
#           (PASS1A / PASS1B preprints, MoTrPAC 2026).
# Inputs  : data/human/705181_target_genes.csv
#             — verified logFC/adj_p for 9 target genes in human skeletal muscle
#               (705181 S3 "Comparisonto_Reitzner" sheet; 13,059-gene complete
#               muscle transcriptome, both modalities × 3 timepoints)
#           data/human/PASS1A_geneset_hits_muscle.csv
#             — Neufer gene set hits from PASS1A 702183 supplements
#           data/human/rat_to_human_gene_mapping.csv
#             — ortholog mapping (rat symbol → human symbol)
# Outputs : figures/fig_human_comparison.png
#           data/human/human_rat_comparison_summary.csv
# Notes   : Raw supplement XLSX files are NOT included in the repository
#           (third-party copyrighted materials). Download from:
#           bioRxiv 10.64898/2026.03.04.705181 (DC4 = S3.xlsx, key file)
#           bioRxiv 10.64898/2026.02.27.702183 (DC5 = PASS1A_S5.xlsx)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
})

# ── Load verified human data ──────────────────────────────────────────────────
target <- read_csv("data/human/705181_target_genes.csv", show_col_types = FALSE)

cat("=== Human MoTrPAC PASS1A — Target Gene Results ===\n")
cat("Source: bioRxiv 10.64898/2026.03.04.705181, S3 Comparisonto_Reitzner sheet\n")
cat("13,059-gene complete muscle transcriptome, n=9 target genes verified\n\n")

print(target, n = Inf)

# ── Key modality × timepoint summary ─────────────────────────────────────────
cat("\n=== ETS genes: selective endurance exercise response (24h) ===\n")
target %>%
  filter(category %in% c("ETS/atp_syn")) %>%
  select(gene, EE_24h_L2FC, EE_24h_adj_p, RE_3.5h_L2FC, RE_3.5h_adj_p) %>%
  print()

cat("\n=== Buffering genes: selective resistance exercise response (3.5h) ===\n")
target %>%
  filter(category %in% c("gsh", "sod", "trx")) %>%
  select(gene, RE_3.5h_L2FC, RE_3.5h_adj_p, EE_24h_L2FC, EE_24h_adj_p) %>%
  print()

cat("\n=== NNT: not regulated by either modality ===\n")
target %>%
  filter(gene == "NNT") %>%
  select(gene, everything()) %>%
  print()

# ── Figure: dot plot of L2FC with significance annotation ─────────────────────
# Key modality × gene comparison
plot_data <- target %>%
  select(gene, category,
         EE_24h_L2FC, EE_24h_adj_p,
         RE_3.5h_L2FC, RE_3.5h_adj_p) %>%
  pivot_longer(
    cols = c(EE_24h_L2FC, RE_3.5h_L2FC),
    names_to = "condition", values_to = "L2FC"
  ) %>%
  mutate(
    adj_p = if_else(condition == "EE_24h_L2FC", EE_24h_adj_p, RE_3.5h_adj_p),
    condition = recode(condition,
      "EE_24h_L2FC"  = "Endurance (24h)",
      "RE_3.5h_L2FC" = "Resistance (3.5h)"),
    significant = adj_p < 0.05,
    gene = factor(gene, levels = rev(sort(unique(gene))))
  )

# Color by functional category
CAT_COL <- c("ETS/atp_syn" = "#2166AC", "gsh" = "#E08214",
             "trx" = "#01665E", "sod" = "#B35806", "nnt" = "#8C510A")

p <- ggplot(plot_data, aes(x = L2FC, y = gene,
                            colour = category, shape = significant)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_point(size = 3, alpha = 0.9, stroke = 1.2) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c("adj_p < 0.05", "n.s."),
                     name = "Significance") +
  scale_colour_manual(values = CAT_COL, name = "Category") +
  facet_wrap(~condition, ncol = 2) +
  labs(
    title = "Human MoTrPAC: target gene response to acute exercise",
    subtitle = "Skeletal muscle, PASS1A (n~20/group). Filled = adj_p < 0.05.",
    x = "log\u2082FC (trained vs. sedentary)", y = NULL,
    caption = "Source: Keshishian et al. bioRxiv 10.64898/2026.03.04.705181, S3"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave("figures/fig_human_comparison.png", p, width = 10, height = 5.5, dpi = 200)
ggsave("figures/fig_human_comparison.pdf", p, width = 10, height = 5.5,
       device = cairo_pdf)
cat("\nSaved figures/fig_human_comparison.png\n")

# ── Summary table ─────────────────────────────────────────────────────────────
target %>%
  write_csv("data/human/human_rat_comparison_summary.csv")
cat("Saved data/human/human_rat_comparison_summary.csv\n")
