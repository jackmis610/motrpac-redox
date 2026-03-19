# =============================================================================
# 06_gene_by_gene.R
# Purpose : Compute pairwise Pearson correlations between all 32 ETS genes and
#           14 buffering genes (cross-tissue 8-week log2FC, males, TRNSCRPT).
#           Produces the 32×14 heatmap showing internal dissociation:
#           GPx4/Txn2 track ETS at r>0.97; NNT trends inversely (r~-0.40).
#           Run with and without BAT; highlights r>0.70 and r<-0.30 pairs.
# Inputs  : data/neufer_da_all.csv   (from 02_query_neufer.R)
# Outputs : figures/fig_gene_by_gene_ets_buffering_heatmap.png/.pdf
#           figures/fig_gene_by_gene_ets_buffering_heatmap_noBat.png/.pdf
#           data/gene_gene_ets_buffering_cor.csv
#           data/gene_gene_ets_buffering_cor_noBat.csv
# =============================================================================
# Fig: Gene-by-gene ETS × buffering correlation heatmap
# Cross-tissue Pearson r of 8-week log2FC (males, TRNSCRPT)
# Output: figures/fig_gene_by_gene_ets_buffering_heatmap.pdf
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ── Gene category definitions ────────────────────────────────
ets_genes <- list(
  CI          = c("Ndufs1","Ndufs2","Ndufs3","Ndufs7","Ndufs8",
                  "Ndufv1","Ndufv2","Ndufa9","Ndufa10","Ndufb8"),
  CIII        = c("Uqcrc1","Uqcrc2","Uqcrfs1","Uqcrb","Cyc1"),
  CIV         = c("Cox4i1","Cox5a","Cox5b","Cox6a1","Cox6b1","Cox7a2"),
  CytC        = c("Cycs"),
  ATPsynthase = c("Atp5f1a","Atp5f1b","Atp5f1c","Atp5f1d","Atp5f1e",
                  "Atp5pb","Atp5mc1","Atp5mc2","Atp5mc3","Atp5po")
)

buf_genes <- list(
  GSH = c("Gpx1","Gpx4","Gsr","Gclc","Gclm","Gss","Glrx2"),
  TRX = c("Txn2","Txnrd2","Prdx3","Prdx5"),
  NNT = c("Nnt"),
  SOD = c("Sod2","Sod1")
)

ets_all  <- unlist(ets_genes,  use.names = FALSE)
buf_all  <- unlist(buf_genes,  use.names = FALSE)

# Lookup: gene → complex label
complex_map <- bind_rows(
  tibble(gene = ets_genes$CI,          complex = "CI"),
  tibble(gene = ets_genes$CIII,        complex = "CIII"),
  tibble(gene = ets_genes$CIV,         complex = "CIV"),
  tibble(gene = ets_genes$CytC,        complex = "Cyt c"),
  tibble(gene = ets_genes$ATPsynthase, complex = "ATP syn")
)

buf_map <- bind_rows(
  tibble(gene = buf_genes$GSH, circuit = "GSH"),
  tibble(gene = buf_genes$TRX, circuit = "TRX"),
  tibble(gene = buf_genes$NNT, circuit = "NNT"),
  tibble(gene = buf_genes$SOD, circuit = "SOD")
)

# ── Load data ────────────────────────────────────────────────
cat("Loading data...\n")
da <- read_csv("data/neufer_da_all.csv", show_col_types = FALSE)

# ── Filter: males, 8w, TRNSCRPT, ETS or buffering ───────────
d8 <- da %>%
  filter(
    assay            == "TRNSCRPT",
    comparison_group == "8w",
    sex              == "male",
    gene_symbol %in% c(ets_all, buf_all)
  )

cat("Rows after filter:", nrow(d8), "\n")

# ── Gene-level mean logFC per tissue (collapse multi-feature) ─
gene_tissue <- d8 %>%
  group_by(tissue, gene_symbol) %>%
  summarise(logFC = mean(logFC, na.rm = TRUE), .groups = "drop")

cat("Tissues present:", n_distinct(gene_tissue$tissue), "\n")
cat("Genes present: ",  n_distinct(gene_tissue$gene_symbol), "\n")

# ── Pivot to tissue × gene matrix ────────────────────────────
mat_wide <- gene_tissue %>%
  pivot_wider(names_from = gene_symbol, values_from = logFC)

tissues <- mat_wide$tissue
mat <- as.matrix(mat_wide[, -1])
rownames(mat) <- tissues

cat("Matrix dimensions:", dim(mat), "\n")

# ── Subset to genes present in data ─────────────────────────
ets_present <- intersect(ets_all, colnames(mat))
buf_present <- intersect(buf_all, colnames(mat))

cat("ETS genes with data:", length(ets_present), "\n")
cat("Buffering genes with data:", length(buf_present), "\n")

mat_ets <- mat[, ets_present, drop = FALSE]
mat_buf <- mat[, buf_present, drop = FALSE]

# Handle tissues with any NA (drop row if NA in either matrix)
complete_rows <- complete.cases(mat_ets) & complete.cases(mat_buf)
mat_ets <- mat_ets[complete_rows, ]
mat_buf <- mat_buf[complete_rows, ]

cat("Tissues used in correlation (complete cases):", sum(complete_rows), "\n")

# ── Compute ETS × buffering correlation matrix ───────────────
# Result: n_ets × n_buf (rows = ETS, cols = buffering)
n_tissues <- nrow(mat_ets)
cor_mat <- matrix(NA_real_,
                  nrow = length(ets_present),
                  ncol = length(buf_present),
                  dimnames = list(ets_present, buf_present))

pval_mat <- matrix(NA_real_,
                   nrow = length(ets_present),
                   ncol = length(buf_present),
                   dimnames = list(ets_present, buf_present))

for (eg in ets_present) {
  for (bg in buf_present) {
    ct <- cor.test(mat_ets[, eg], mat_buf[, bg],
                   method = "pearson", use = "complete.obs")
    cor_mat[eg, bg]  <- ct$estimate
    pval_mat[eg, bg] <- ct$p.value
  }
}

# ── Report top and bottom pairs ──────────────────────────────
cor_long <- as.data.frame(as.table(cor_mat)) %>%
  rename(ets_gene = Var1, buf_gene = Var2, r = Freq) %>%
  arrange(desc(r))

cat("\n=== TOP 10 POSITIVE ETS-BUFFERING GENE PAIRS ===\n")
print(head(cor_long, 10))

cat("\n=== TOP 10 NEGATIVE ETS-BUFFERING GENE PAIRS ===\n")
print(tail(cor_long, 10) %>% arrange(r))

write_csv(cor_long %>% arrange(desc(r)),
          "data/gene_gene_ets_buffering_cor.csv")
cat("\nSaved data/gene_gene_ets_buffering_cor.csv\n")

# ── Hierarchical clustering ───────────────────────────────────
clust_ets <- hclust(dist(cor_mat),     method = "ward.D2")
clust_buf <- hclust(dist(t(cor_mat)),  method = "ward.D2")

# Reorder
cor_mat_ord <- cor_mat[clust_ets$order, clust_buf$order]

# ── Annotations ──────────────────────────────────────────────
# ETS row annotation: complex
ets_ord_names <- rownames(cor_mat_ord)
complex_labels <- complex_map$complex[match(ets_ord_names, complex_map$gene)]

complex_colors <- c(
  "CI"      = "#2166AC",
  "CIII"    = "#4DAC26",
  "CIV"     = "#D01C8B",
  "Cyt c"   = "#F1A340",
  "ATP syn" = "#762A83"
)

row_ann <- rowAnnotation(
  Complex = complex_labels,
  col = list(Complex = complex_colors),
  annotation_name_gp = gpar(fontsize = 9),
  annotation_legend_param = list(
    Complex = list(title = "ETS complex", title_gp = gpar(fontsize = 9),
                   labels_gp = gpar(fontsize = 8))
  )
)

# Buffering col annotation: circuit
buf_ord_names <- colnames(cor_mat_ord)
circuit_labels <- buf_map$circuit[match(buf_ord_names, buf_map$gene)]

circuit_colors <- c(
  "GSH" = "#E08214",
  "TRX" = "#01665E",
  "NNT" = "#8C510A",
  "SOD" = "#C7EAE5"
)

col_ann <- HeatmapAnnotation(
  Circuit = circuit_labels,
  col = list(Circuit = circuit_colors),
  annotation_name_gp = gpar(fontsize = 9),
  annotation_legend_param = list(
    Circuit = list(title = "Buffering circuit", title_gp = gpar(fontsize = 9),
                   labels_gp = gpar(fontsize = 8))
  )
)

# ── Cell function: highlight r > 0.7 or r < -0.3 ───────────
cell_fn <- function(j, i, x, y, width, height, fill) {
  v <- cor_mat_ord[i, j]
  grid.rect(x, y, width, height,
            gp = gpar(col = NA, fill = fill))
  if (!is.na(v) && (v > 0.7 || v < -0.3)) {
    # Bold border
    border_col <- if (v > 0.7) "black" else "#8B0000"
    lwd        <- if (v > 0.7) 1.8 else 1.4
    grid.rect(x, y, width, height,
              gp = gpar(col = border_col, fill = NA, lwd = lwd))
  }
}

# ── Correlation color scale ───────────────────────────────────
col_fun <- colorRamp2(
  c(-1, -0.5, 0, 0.5, 1),
  c("#D01C8B", "#F1B6DA", "white", "#B8E186", "#4DAC26")
)

# ── Build heatmap ─────────────────────────────────────────────
ht <- Heatmap(
  cor_mat_ord,
  name                  = "Pearson r",
  col                   = col_fun,
  cluster_rows          = FALSE,   # already ordered
  cluster_columns       = FALSE,
  row_dend_reorder      = FALSE,
  column_dend_reorder   = FALSE,
  left_annotation       = row_ann,
  top_annotation        = col_ann,
  cell_fun              = cell_fn,
  row_names_gp          = gpar(fontsize = 7.5),
  column_names_gp       = gpar(fontsize = 7.5),
  column_names_rot      = 45,
  row_title             = "ETS genes (clustered)",
  column_title          = "Redox buffering genes (clustered)",
  row_title_gp          = gpar(fontsize = 10, fontface = "bold"),
  column_title_gp       = gpar(fontsize = 10, fontface = "bold"),
  heatmap_legend_param  = list(
    title     = "Pearson r",
    title_gp  = gpar(fontsize = 9),
    labels_gp = gpar(fontsize = 8),
    at        = c(-1, -0.5, 0, 0.5, 1),
    direction = "vertical"
  ),
  width  = unit(ncol(cor_mat_ord) * 0.55, "cm"),
  height = unit(nrow(cor_mat_ord) * 0.55, "cm")
)

# ── Save ─────────────────────────────────────────────────────
pdf_file <- "figures/fig_gene_by_gene_ets_buffering_heatmap.pdf"
png_file <- "figures/fig_gene_by_gene_ets_buffering_heatmap.png"

pdf(pdf_file, width = 14, height = 12)
draw(ht,
     column_title      = "Cross-tissue Pearson r: ETS vs. Redox Buffering (Males, 8w TRNSCRPT)",
     column_title_gp   = gpar(fontsize = 12, fontface = "bold"),
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     padding = unit(c(5, 5, 5, 5), "mm"))

# Legend note for highlighted borders
grid.text(
  "Bold black border: r > 0.70 | Dark red border: r < -0.30",
  x = unit(0.5, "npc"), y = unit(0.01, "npc"),
  gp = gpar(fontsize = 7.5, col = "grey40", fontface = "italic"),
  just = "center"
)
dev.off()

png(png_file, width = 2800, height = 2400, res = 200)
draw(ht,
     column_title      = "Cross-tissue Pearson r: ETS vs. Redox Buffering (Males, 8w TRNSCRPT)",
     column_title_gp   = gpar(fontsize = 12, fontface = "bold"),
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     padding = unit(c(5, 5, 5, 5), "mm"))
grid.text(
  "Bold black border: r > 0.70 | Dark red border: r < -0.30",
  x = unit(0.5, "npc"), y = unit(0.01, "npc"),
  gp = gpar(fontsize = 7.5, col = "grey40", fontface = "italic"),
  just = "center"
)
dev.off()

cat("\nSaved", pdf_file, "\n")
cat("Saved", png_file, "\n")

# ── Summary statistics ────────────────────────────────────────
cat("\n=== SUMMARY ===\n")
cat("ETS genes in matrix:", nrow(cor_mat), "\n")
cat("Buffering genes in matrix:", ncol(cor_mat), "\n")
cat("Tissues used:", n_tissues, "\n")
cat("Pairs r > 0.70:", sum(cor_mat > 0.70, na.rm = TRUE), "\n")
cat("Pairs r < -0.30:", sum(cor_mat < -0.30, na.rm = TRUE), "\n")
cat("Overall mean r:", round(mean(cor_mat, na.rm = TRUE), 3), "\n")
cat("Overall median r:", round(median(cor_mat, na.rm = TRUE), 3), "\n")

cat("\n=== HIGHEST SINGLE PAIR ===\n")
top1 <- cor_long[1, ]
cat(sprintf("  %s ~ %s: r = %.3f\n", top1$ets_gene, top1$buf_gene, top1$r))

cat("\n=== LOWEST SINGLE PAIR ===\n")
bot1 <- cor_long[nrow(cor_long), ]
cat(sprintf("  %s ~ %s: r = %.3f\n", bot1$ets_gene, bot1$buf_gene, bot1$r))

cat("\nDone.\n")
