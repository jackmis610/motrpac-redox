# ============================================================
# Sensitivity analysis: all three checks excluding BAT
#   1. Main ETS vs buffering regression (both sexes)
#   2. Permutation test (1,000 permutations) — R² specificity
#   3. Gene-by-gene Pearson correlation matrix (32 × 14)
#
# All outputs saved with _noBat suffix
# ============================================================

suppressPackageStartupMessages({
  library(MotrpacRatTraining6moData)
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

set.seed(2026)
EXCLUDE_TISSUE <- "BAT"

cat("=== Sensitivity analysis: excluding", EXCLUDE_TISSUE, "===\n\n")

# ── Shared gene categories ────────────────────────────────────
ets_cats <- c("ets_CI","ets_CIII","ets_CIV","cyt_c","atp_synthase")
buf_cats <- c("gsh","trx","nnt","sod")
SEX_COL  <- c(Male = "#2166AC", Female = "#D6604D")

# ── Load curated DA data ─────────────────────────────────────
cat("[Data] Loading neufer_da_all.csv...\n")
da <- read_csv("data/neufer_da_all.csv", show_col_types = FALSE)

# ── Observed values WITH BAT (for comparison) ────────────────
obs_with_bat <- list(
  male_slope = 0.6507, male_r2 = 0.7163, male_p = 1.78e-03,
  fem_slope  = 0.1497, fem_r2  = 0.3331, fem_p  = 1.38e-11
)

# ============================================================
# SECTION 1: Main regression without BAT
# ============================================================
cat("\n============================================================\n")
cat("SECTION 1: Main ETS vs buffering regression (no BAT)\n")
cat("============================================================\n")

make_indices_noBAT <- function(sex_filter) {
  da %>%
    filter(assay == "TRNSCRPT",
           comparison_group == "8w",
           sex == sex_filter,
           tissue != EXCLUDE_TISSUE,
           framework_category %in% c(ets_cats, buf_cats)) %>%
    mutate(gene_class = if_else(framework_category %in% ets_cats,
                                "ets", "buffering")) %>%
    group_by(tissue, gene_symbol, gene_class) %>%
    summarise(gene_logFC = mean(logFC, na.rm = TRUE), .groups = "drop") %>%
    group_by(tissue, gene_class) %>%
    summarise(index_mean = mean(gene_logFC, na.rm = TRUE),
              index_se   = sd(gene_logFC, na.rm = TRUE) / sqrt(n()),
              .groups = "drop") %>%
    pivot_wider(names_from = gene_class,
                values_from = c(index_mean, index_se)) %>%
    rename(ets_mean = index_mean_ets,
           buf_mean = index_mean_buffering,
           ets_se   = index_se_ets,
           buf_se   = index_se_buffering) %>%
    filter(!is.na(ets_mean), !is.na(buf_mean))
}

idx_m <- make_indices_noBAT("male")
idx_f <- make_indices_noBAT("female")
cat("Tissues used (no BAT):", nrow(idx_m), "\n")

fit_m <- lm(buf_mean ~ ets_mean, data = idx_m)
fit_f <- lm(buf_mean ~ ets_mean, data = idx_f)

get_stats <- function(fit) {
  s  <- summary(fit)
  b  <- coef(fit)[2]
  se <- s$coefficients["ets_mean","Std. Error"]
  n  <- nrow(fit$model)
  list(slope = b, r2 = s$r.squared, intercept = coef(fit)[1],
       p_vs1 = pt((b - 1)/se, df = n - 2), n = n, se = se)
}

sm <- get_stats(fit_m)
sf <- get_stats(fit_f)

results_reg <- tibble(
  sex         = c("Male","Female"),
  n_tissues   = c(sm$n, sf$n),
  slope       = c(sm$slope, sf$slope),
  r_squared   = c(sm$r2, sf$r2),
  p_slope_lt1 = c(sm$p_vs1, sf$p_vs1),
  slope_WITH_bat = c(obs_with_bat$male_slope, obs_with_bat$fem_slope),
  r2_WITH_bat    = c(obs_with_bat$male_r2,    obs_with_bat$fem_r2),
  p_WITH_bat     = c(obs_with_bat$male_p,     obs_with_bat$fem_p)
)

cat("\n=== REGRESSION RESULTS (no BAT) ===\n")
results_reg %>%
  mutate(across(c(slope, r_squared, slope_WITH_bat, r2_WITH_bat), ~round(.x, 4)),
         across(c(p_slope_lt1, p_WITH_bat), ~signif(.x, 3))) %>%
  print()

write_csv(results_reg, "data/regression_noBat.csv")
cat("Saved data/regression_noBat.csv\n")

# ── Scatter plot ──────────────────────────────────────────────
make_scatter <- function(idx, stats, sex_label, col) {
  xs <- seq(min(idx$ets_mean), max(idx$ets_mean), length.out = 50)
  line_df <- tibble(x = xs, y = stats$intercept + stats$slope * xs)
  ann <- sprintf("\u03b2 = %.3f (was %.3f)\nR\u00b2 = %.3f (was %.3f)\np(\u03b2<1) = %s",
                 stats$slope, if(sex_label=="Male") obs_with_bat$male_slope else obs_with_bat$fem_slope,
                 stats$r2,    if(sex_label=="Male") obs_with_bat$male_r2    else obs_with_bat$fem_r2,
                 if(stats$p_vs1 < 0.001) "<0.001" else sprintf("%.3f", stats$p_vs1))
  ggplot(idx, aes(ets_mean, buf_mean)) +
    geom_abline(slope=1, intercept=0, linetype="dashed",
                colour="grey60", linewidth=0.5) +
    geom_hline(yintercept=0, colour="grey85", linewidth=0.3) +
    geom_vline(xintercept=0, colour="grey85", linewidth=0.3) +
    geom_line(data=line_df, aes(x=x,y=y), colour=col,
              linewidth=0.9, inherit.aes=FALSE) +
    geom_errorbar(aes(ymin=buf_mean-buf_se, ymax=buf_mean+buf_se),
                  width=0, alpha=0.35, colour=col, linewidth=0.4) +
    geom_errorbar(aes(xmin=ets_mean-ets_se, xmax=ets_mean+ets_se),
                  orientation="y", width=0, alpha=0.35,
                  colour=col, linewidth=0.4) +
    geom_point(colour=col, size=2.5, alpha=0.9) +
    geom_text_repel(aes(label=tissue), size=2.2, colour=col,
                    max.overlaps=16, segment.colour="grey60",
                    segment.size=0.3) +
    annotate("text", x=Inf, y=-Inf, hjust=1.05, vjust=-0.4,
             label=ann, size=2.9, fontface="italic", colour=col) +
    labs(title=sprintf("%s — ETS vs buffering (BAT excluded)", sex_label),
         x="ETS index (mean log\u2082FC)",
         y="Buffering index (mean log\u2082FC)") +
    theme_classic(base_size=10) +
    theme(plot.title=element_text(face="bold", size=10, hjust=0.5))
}

p_reg <- make_scatter(idx_m, sm, "Male",   SEX_COL["Male"]) +
         make_scatter(idx_f, sf, "Female", SEX_COL["Female"])

p_reg_full <- p_reg +
  plot_annotation(
    title    = "ETS vs redox buffering scaling — BAT excluded (8w, TRNSCRPT)",
    subtitle = "Dashed line = 1:1. Values in parentheses = full model including BAT.",
    theme = theme(
      plot.title    = element_text(face="bold", size=12, hjust=0.5),
      plot.subtitle = element_text(colour="grey40", size=8.5, hjust=0.5)
    )
  )

ggsave("figures/fig_regression_noBat.pdf",
       p_reg_full, width=12, height=6, device=cairo_pdf)
ggsave("figures/fig_regression_noBat.png",
       p_reg_full, width=12, height=6, dpi=200)
cat("Saved figures/fig_regression_noBat.pdf/.png\n")

# ============================================================
# SECTION 2: Permutation test without BAT
# ============================================================
cat("\n============================================================\n")
cat("SECTION 2: Permutation test (no BAT)\n")
cat("============================================================\n")

N_PERM <- 1000L
N_ETS  <- 32L
N_BUF  <- 14L
MIN_T  <- 11L   # one less than before (17 tissues now)

cat("[Data] Loading full transcriptome for permutation...\n")
data(FEATURE_TO_GENE)
f2g_min <- FEATURE_TO_GENE %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  select(feature_ID, gene_symbol) %>%
  distinct()

all_items   <- data(package = "MotrpacRatTraining6moData")$results[,"Item"]
trnscrpt_da <- all_items[grepl("^TRNSCRPT_.*_DA$", all_items)]

all_logfc <- map_dfr(trnscrpt_da, function(tn) {
  e <- new.env()
  tryCatch({
    data(list=tn, package="MotrpacRatTraining6moData", envir=e)
    get(tn, envir=e) %>%
      filter(sex=="male", comparison_group=="8w",
             tissue != EXCLUDE_TISSUE) %>%
      select(feature_ID, tissue, logFC) %>%
      inner_join(f2g_min, by="feature_ID") %>%
      group_by(tissue, gene_symbol) %>%
      summarise(logFC = mean(logFC, na.rm=TRUE), .groups="drop")
  }, error = function(e) tibble())
})

cat("Tissues after BAT exclusion:", n_distinct(all_logfc$tissue), "\n")

gene_cov <- all_logfc %>%
  group_by(gene_symbol) %>%
  summarise(n_t = n_distinct(tissue), .groups="drop")

universe <- gene_cov %>% filter(n_t >= MIN_T) %>% pull(gene_symbol)
cat("Gene universe (≥", MIN_T, "tissues):", length(universe), "\n")

mat_wide <- all_logfc %>%
  filter(gene_symbol %in% universe) %>%
  group_by(tissue, gene_symbol) %>%
  summarise(logFC=mean(logFC, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from=gene_symbol, values_from=logFC, values_fill=NA_real_)

mat <- as.matrix(mat_wide[,-1])
rownames(mat) <- mat_wide$tissue
cat("Permutation matrix:", nrow(mat), "tissues x", ncol(mat), "genes\n")

# Observed slope and R² without BAT (males)
source("scripts/neufer_geneset.R", local=TRUE)
ets_obs <- intersect(unlist(neufer_geneset[ets_cats]), colnames(mat))
buf_obs <- intersect(unlist(neufer_geneset[buf_cats]), colnames(mat))

ets_m_obs <- rowMeans(mat[, ets_obs, drop=FALSE], na.rm=TRUE)
buf_m_obs <- rowMeans(mat[, buf_obs, drop=FALSE], na.rm=TRUE)
obs_fit   <- lm(buf_m_obs ~ ets_m_obs)
obs_slope_nb <- coef(obs_fit)[2]
obs_r2_nb    <- summary(obs_fit)$r.squared

cat(sprintf("\nObserved (no BAT): slope=%.4f  R²=%.4f\n", obs_slope_nb, obs_r2_nb))

# Complete genes (no NA)
complete_genes <- colnames(mat)[colSums(is.na(mat)) == 0]
cat("Genes with no NA:", length(complete_genes), "\n")

if (length(complete_genes) < (N_ETS + N_BUF) * 3) {
  complete_genes <- gene_cov %>%
    filter(n_t >= nrow(mat)) %>%
    pull(gene_symbol) %>%
    intersect(colnames(mat))
  cat("Relaxed: complete genes:", length(complete_genes), "\n")
}

perm_mat <- mat[, complete_genes, drop=FALSE]

null_stats <- matrix(NA_real_, nrow=N_PERM, ncol=2,
                     dimnames=list(NULL, c("slope","r2")))
cat(sprintf("Running %d permutations...\n", N_PERM))
for (i in seq_len(N_PERM)) {
  if (i %% 200 == 0) cat("  perm", i, "\n")
  sel <- sample(complete_genes, N_ETS + N_BUF, replace=FALSE)
  em <- rowMeans(perm_mat[, sel[1:N_ETS],            drop=FALSE], na.rm=TRUE)
  bm <- rowMeans(perm_mat[, sel[(N_ETS+1):(N_ETS+N_BUF)], drop=FALSE], na.rm=TRUE)
  ok <- !is.na(em) & !is.na(bm)
  if (sum(ok) < 4) next
  fit <- lm(bm[ok] ~ em[ok])
  null_stats[i,] <- c(coef(fit)[2], summary(fit)$r.squared)
}

null_stats <- null_stats[complete.cases(null_stats),]
null_slopes_nb <- null_stats[,"slope"]
null_r2_nb     <- null_stats[,"r2"]

p_slope_nb <- mean(null_slopes_nb <= obs_slope_nb)
p_r2_nb    <- mean(null_r2_nb    >= obs_r2_nb)

cat(sprintf("\n=== PERMUTATION RESULTS (no BAT) ===\n"))
cat(sprintf("  Observed slope: %.4f  (with BAT: %.4f)\n", obs_slope_nb, 0.6507))
cat(sprintf("  Observed R²:    %.4f  (with BAT: %.4f)\n", obs_r2_nb,    0.7163))
cat(sprintf("  Null mean R²:   %.4f\n", mean(null_r2_nb)))
cat(sprintf("  p(slope≤obs):   %.4f  (with BAT: 0.8100)\n", p_slope_nb))
cat(sprintf("  p(R²≥obs):      %.4f  (with BAT: 0.0200) [specificity signal]\n", p_r2_nb))

write_csv(as_tibble(null_stats),
          "data/permutation_null_slopes_noBat.csv")
cat("Saved data/permutation_null_slopes_noBat.csv\n")

# ── Plot permutation ──────────────────────────────────────────
fmt_p <- function(p) if(p<0.001) "p < 0.001" else sprintf("p = %.3f", p)
null_df <- tibble(slope=null_slopes_nb, r2=null_r2_nb)

pP_slope <- ggplot(null_df, aes(x=slope)) +
  geom_histogram(aes(y=after_stat(density)), bins=60,
                 fill="grey75", colour="white", linewidth=0.2) +
  geom_density(colour="grey40", linewidth=0.7) +
  geom_vline(xintercept=obs_slope_nb, colour="#D01C8B",
             linewidth=1.3) +
  annotate("rect", xmin=-Inf, xmax=obs_slope_nb,
           ymin=-Inf, ymax=Inf, fill="#D01C8B", alpha=0.07) +
  annotate("text", x=obs_slope_nb+0.04, y=Inf,
           label=sprintf("Observed \u03b2=%.3f\n%s", obs_slope_nb, fmt_p(p_slope_nb)),
           hjust=0, vjust=1.3, colour="#D01C8B", size=3.2, fontface="bold") +
  labs(title="A. Slope (\u03b2) — no BAT",
       x="Regression slope (\u03b2)", y="Density") +
  theme_classic(base_size=10) +
  theme(plot.title=element_text(face="bold", size=10))

pP_r2 <- ggplot(null_df, aes(x=r2)) +
  geom_histogram(aes(y=after_stat(density)), bins=60,
                 fill="grey75", colour="white", linewidth=0.2) +
  geom_density(colour="grey40", linewidth=0.7) +
  geom_vline(xintercept=obs_r2_nb, colour="#2166AC",
             linewidth=1.3) +
  annotate("rect", xmin=obs_r2_nb, xmax=Inf,
           ymin=-Inf, ymax=Inf, fill="#2166AC", alpha=0.09) +
  annotate("text", x=obs_r2_nb-0.02, y=Inf,
           label=sprintf("Observed R\u00b2=%.3f\n%s", obs_r2_nb, fmt_p(p_r2_nb)),
           hjust=1, vjust=1.3, colour="#2166AC", size=3.2, fontface="bold") +
  labs(title="B. R\u00b2 (coupling) — no BAT",
       x="Variance explained (R\u00b2)", y="Density") +
  theme_classic(base_size=10) +
  theme(plot.title=element_text(face="bold", size=10))

pP <- pP_slope + pP_r2 +
  plot_annotation(
    title    = "Permutation test (BAT excluded) — 1,000 random gene sets, males, 8w",
    subtitle = sprintf("Universe: %s genes, 17 tissues",
                       format(length(complete_genes), big.mark=",")),
    theme = theme(
      plot.title    = element_text(face="bold", size=12, hjust=0.5),
      plot.subtitle = element_text(colour="grey40", size=8.5, hjust=0.5)
    )
  )

ggsave("figures/fig_permutation_null_noBat.pdf",
       pP, width=12, height=5.5, device=cairo_pdf)
ggsave("figures/fig_permutation_null_noBat.png",
       pP, width=12, height=5.5, dpi=200)
cat("Saved figures/fig_permutation_null_noBat.pdf/.png\n")

# ============================================================
# SECTION 3: Gene-by-gene correlation matrix without BAT
# ============================================================
cat("\n============================================================\n")
cat("SECTION 3: Gene-by-gene correlation matrix (no BAT)\n")
cat("============================================================\n")

ets_gene_list <- list(
  CI          = c("Ndufs1","Ndufs2","Ndufs3","Ndufs7","Ndufs8",
                  "Ndufv1","Ndufv2","Ndufa9","Ndufa10","Ndufb8"),
  CIII        = c("Uqcrc1","Uqcrc2","Uqcrfs1","Uqcrb","Cyc1"),
  CIV         = c("Cox4i1","Cox5a","Cox5b","Cox6a1","Cox6b1","Cox7a2"),
  CytC        = c("Cycs"),
  ATPsynthase = c("Atp5f1a","Atp5f1b","Atp5f1c","Atp5f1d","Atp5f1e",
                  "Atp5pb","Atp5mc1","Atp5mc2","Atp5mc3","Atp5po")
)
buf_gene_list <- list(
  GSH = c("Gpx1","Gpx4","Gsr","Gclc","Gclm","Gss","Glrx2"),
  TRX = c("Txn2","Txnrd2","Prdx3","Prdx5"),
  NNT = c("Nnt"),
  SOD = c("Sod2","Sod1")
)

ets_all  <- unlist(ets_gene_list, use.names=FALSE)
buf_all  <- unlist(buf_gene_list, use.names=FALSE)
cplx_map <- bind_rows(
  tibble(gene=ets_gene_list$CI,          complex="CI"),
  tibble(gene=ets_gene_list$CIII,        complex="CIII"),
  tibble(gene=ets_gene_list$CIV,         complex="CIV"),
  tibble(gene=ets_gene_list$CytC,        complex="Cyt c"),
  tibble(gene=ets_gene_list$ATPsynthase, complex="ATP syn")
)
buf_map <- bind_rows(
  tibble(gene=buf_gene_list$GSH, circuit="GSH"),
  tibble(gene=buf_gene_list$TRX, circuit="TRX"),
  tibble(gene=buf_gene_list$NNT, circuit="NNT"),
  tibble(gene=buf_gene_list$SOD, circuit="SOD")
)

# Build gene × tissue matrix from neufer_da_all (no BAT)
gene_tissue_nb <- da %>%
  filter(assay=="TRNSCRPT", comparison_group=="8w", sex=="male",
         tissue != EXCLUDE_TISSUE,
         gene_symbol %in% c(ets_all, buf_all)) %>%
  group_by(tissue, gene_symbol) %>%
  summarise(logFC=mean(logFC, na.rm=TRUE), .groups="drop")

mat_nb_wide <- gene_tissue_nb %>%
  pivot_wider(names_from=gene_symbol, values_from=logFC)

mat_nb <- as.matrix(mat_nb_wide[,-1])
rownames(mat_nb) <- mat_nb_wide$tissue

ets_in <- intersect(ets_all, colnames(mat_nb))
buf_in <- intersect(buf_all, colnames(mat_nb))

mat_ets_nb <- mat_nb[, ets_in, drop=FALSE]
mat_buf_nb <- mat_nb[, buf_in, drop=FALSE]
ok_rows    <- complete.cases(mat_ets_nb) & complete.cases(mat_buf_nb)
mat_ets_nb <- mat_ets_nb[ok_rows,]
mat_buf_nb <- mat_buf_nb[ok_rows,]

cat("Tissues in correlation (no BAT):", sum(ok_rows), "\n")

# Compute correlations
cor_nb <- matrix(NA_real_, nrow=length(ets_in), ncol=length(buf_in),
                 dimnames=list(ets_in, buf_in))
for (eg in ets_in) for (bg in buf_in) {
  ct <- cor.test(mat_ets_nb[,eg], mat_buf_nb[,bg],
                 method="pearson", use="complete.obs")
  cor_nb[eg, bg] <- ct$estimate
}

# Also load original for comparison
cor_orig <- read_csv("data/gene_gene_ets_buffering_cor.csv",
                     show_col_types=FALSE) %>%
  rename(r_with_bat = r)

cor_nb_long <- as.data.frame(as.table(cor_nb)) %>%
  rename(ets_gene=Var1, buf_gene=Var2, r=Freq) %>%
  as_tibble() %>%
  left_join(cor_orig, by=c("ets_gene","buf_gene")) %>%
  mutate(delta_r = r - r_with_bat) %>%
  arrange(desc(r))

write_csv(cor_nb_long, "data/gene_gene_ets_buffering_cor_noBat.csv")
cat("Saved data/gene_gene_ets_buffering_cor_noBat.csv\n")

cat("\n=== TOP 5 POSITIVE PAIRS (no BAT) ===\n")
cor_nb_long %>% head(5) %>%
  select(ets_gene, buf_gene, r, r_with_bat, delta_r) %>%
  mutate(across(where(is.numeric), ~round(.x, 4))) %>% print()

cat("\n=== TOP 5 NEGATIVE PAIRS (no BAT) ===\n")
cor_nb_long %>% arrange(r) %>% head(5) %>%
  select(ets_gene, buf_gene, r, r_with_bat, delta_r) %>%
  mutate(across(where(is.numeric), ~round(.x, 4))) %>% print()

cat("\n=== NNT pairs specifically (no BAT) ===\n")
cor_nb_long %>%
  filter(buf_gene == "Nnt") %>%
  arrange(r) %>%
  select(ets_gene, buf_gene, r, r_with_bat, delta_r) %>%
  mutate(across(where(is.numeric), ~round(.x, 4))) %>%
  print(n=20)

cat(sprintf("\nNNT pairs still negative: %d / %d\n",
            sum(cor_nb_long$buf_gene=="Nnt" & cor_nb_long$r < 0, na.rm=TRUE),
            sum(cor_nb_long$buf_gene=="Nnt", na.rm=TRUE)))
cat(sprintf("NNT pairs r < -0.3: %d / %d\n",
            sum(cor_nb_long$buf_gene=="Nnt" & cor_nb_long$r < -0.3, na.rm=TRUE),
            sum(cor_nb_long$buf_gene=="Nnt", na.rm=TRUE)))

# Summary stats
cat(sprintf("\nPairs r > 0.70: %d (was %d with BAT)\n",
            sum(cor_nb > 0.70, na.rm=TRUE), 181))
cat(sprintf("Pairs r < -0.30: %d (was %d with BAT)\n",
            sum(cor_nb < -0.30, na.rm=TRUE), 60))
cat(sprintf("Overall mean r: %.3f (was 0.435 with BAT)\n",
            mean(cor_nb, na.rm=TRUE)))

# ── Heatmap ───────────────────────────────────────────────────
clust_e <- hclust(dist(cor_nb),    method="ward.D2")
clust_b <- hclust(dist(t(cor_nb)), method="ward.D2")
cor_ord  <- cor_nb[clust_e$order, clust_b$order]

e_ord <- rownames(cor_ord)
b_ord <- colnames(cor_ord)

cplx_colors   <- c("CI"="chocolate","CIII"="#4DAC26","CIV"="#D01C8B",
                   "Cyt c"="#F1A340","ATP syn"="#762A83")
circuit_colors <- c("GSH"="#E08214","TRX"="#01665E","NNT"="#8C510A","SOD"="#C7EAE5")

row_ann_nb <- rowAnnotation(
  Complex = cplx_map$complex[match(e_ord, cplx_map$gene)],
  col = list(Complex = cplx_colors),
  annotation_name_gp = gpar(fontsize=9)
)
col_ann_nb <- HeatmapAnnotation(
  Circuit = buf_map$circuit[match(b_ord, buf_map$gene)],
  col = list(Circuit = circuit_colors),
  annotation_name_gp = gpar(fontsize=9)
)

col_fun <- colorRamp2(c(-1,-0.5,0,0.5,1),
                      c("#D01C8B","#F1B6DA","white","#B8E186","#4DAC26"))

cell_fn_nb <- function(j,i,x,y,width,height,fill) {
  v <- cor_ord[i,j]
  grid.rect(x,y,width,height,gp=gpar(col=NA,fill=fill))
  if (!is.na(v) && (v>0.7 || v<(-0.3))) {
    bc  <- if(v>0.7) "black" else "#8B0000"
    lwd <- if(v>0.7) 1.8 else 1.4
    grid.rect(x,y,width,height,gp=gpar(col=bc,fill=NA,lwd=lwd))
  }
}

ht_nb <- Heatmap(
  cor_ord, name="Pearson r", col=col_fun,
  cluster_rows=FALSE, cluster_columns=FALSE,
  left_annotation=row_ann_nb, top_annotation=col_ann_nb,
  cell_fun=cell_fn_nb,
  row_names_gp    = gpar(fontsize=7.5),
  column_names_gp = gpar(fontsize=7.5),
  column_names_rot = 45,
  row_title    = "ETS genes (clustered)",
  column_title = "Redox buffering genes (clustered)",
  row_title_gp    = gpar(fontsize=10, fontface="bold"),
  column_title_gp = gpar(fontsize=10, fontface="bold"),
  heatmap_legend_param = list(
    title="Pearson r", title_gp=gpar(fontsize=9),
    labels_gp=gpar(fontsize=8), at=c(-1,-0.5,0,0.5,1)
  ),
  width  = unit(ncol(cor_ord)*0.55, "cm"),
  height = unit(nrow(cor_ord)*0.55, "cm")
)

pdf("figures/fig_gene_by_gene_ets_buffering_heatmap_noBat.pdf",
    width=14, height=12)
draw(ht_nb,
     column_title    = "Cross-tissue Pearson r: ETS vs Buffering (Males, 8w, BAT excluded)",
     column_title_gp = gpar(fontsize=12, fontface="bold"),
     padding         = unit(c(5,5,5,5),"mm"))
grid.text("Bold black: r>0.70 | Dark red: r<-0.30",
          x=unit(0.5,"npc"), y=unit(0.01,"npc"),
          gp=gpar(fontsize=7.5, col="grey40", fontface="italic"),
          just="center")
dev.off()

png("figures/fig_gene_by_gene_ets_buffering_heatmap_noBat.png",
    width=2800, height=2400, res=200)
draw(ht_nb,
     column_title    = "Cross-tissue Pearson r: ETS vs Buffering (Males, 8w, BAT excluded)",
     column_title_gp = gpar(fontsize=12, fontface="bold"),
     padding         = unit(c(5,5,5,5),"mm"))
grid.text("Bold black: r>0.70 | Dark red: r<-0.30",
          x=unit(0.5,"npc"), y=unit(0.01,"npc"),
          gp=gpar(fontsize=7.5, col="grey40", fontface="italic"),
          just="center")
dev.off()

cat("\nSaved figures/fig_gene_by_gene_ets_buffering_heatmap_noBat.pdf/.png\n")

# ============================================================
# SUMMARY TABLE
# ============================================================
cat("\n============================================================\n")
cat("SUMMARY: WITH BAT vs WITHOUT BAT\n")
cat("============================================================\n")

cat(sprintf("\n%-40s %10s %10s\n", "Metric", "With BAT", "No BAT"))
cat(strrep("-", 62), "\n")
cat(sprintf("%-40s %10s %10s\n", "Male slope",
            sprintf("%.3f", obs_with_bat$male_slope),
            sprintf("%.3f", sm$slope)))
cat(sprintf("%-40s %10s %10s\n", "Male R²",
            sprintf("%.3f", obs_with_bat$male_r2),
            sprintf("%.3f", sm$r2)))
cat(sprintf("%-40s %10s %10s\n", "Male p(β<1)",
            sprintf("%.3e", obs_with_bat$male_p),
            sprintf("%.3e", sm$p_vs1)))
cat(sprintf("%-40s %10s %10s\n", "Female slope",
            sprintf("%.3f", obs_with_bat$fem_slope),
            sprintf("%.3f", sf$slope)))
cat(sprintf("%-40s %10s %10s\n", "Female R²",
            sprintf("%.3f", obs_with_bat$fem_r2),
            sprintf("%.3f", sf$r2)))
cat(sprintf("%-40s %10s %10s\n", "Female p(β<1)",
            sprintf("%.3e", obs_with_bat$fem_p),
            sprintf("%.3e", sf$p_vs1)))
cat(sprintf("%-40s %10s %10s\n", "Perm: p(slope≤obs)",
            "0.810", sprintf("%.3f", p_slope_nb)))
cat(sprintf("%-40s %10s %10s\n", "Perm: p(R²≥obs) [specificity]",
            "0.020", sprintf("%.3f", p_r2_nb)))
cat(sprintf("%-40s %10s %10s\n", "Gene-gene: r>0.70 pairs",
            "181", sum(cor_nb > 0.70, na.rm=TRUE)))
cat(sprintf("%-40s %10s %10s\n", "Gene-gene: NNT negative pairs",
            "32/32", sprintf("%d/32", sum(cor_nb_long$buf_gene=="Nnt" & cor_nb_long$r<0))))

cat("\nDone.\n")
