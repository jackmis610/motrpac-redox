# ============================================================
# 06: Targeted follow-ups
#   (1) DRP1 phosphoproteomics — Dnm1l site-level regulation
#   (2) Lactate axis tissue coordination — producer/consumer map
#   (3) Succinate-SDH in lung — metabolite vs. enzyme trajectory
# All figures saved with prefix fig_followup_
# ============================================================

suppressPackageStartupMessages({
  library(MotrpacRatTraining6moData)
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
})

dir.create("figures", showWarnings = FALSE)
dir.create("data",    showWarnings = FALSE)

SEX_COL   <- c(Male = "#2166AC", Female = "#D6604D")
SEX_SHAPE <- c(Male = 16L, Female = 17L)
DIVPAL    <- "RdBu"
TP_LEVELS <- c("1w","2w","4w","8w")

TISSUE_ORDER <- c("SKM-GN","SKM-VL","HEART","BAT","WAT-SC",
                  "LIVER","KIDNEY","LUNG","ADRNL","BLOOD",
                  "CORTEX","HIPPOC","HYPOTH","VENACV",
                  "COLON","SMLINT","SPLEEN","OVARY","TESTES")

theme_pub <- theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text       = element_text(face = "bold", size = 9),
        panel.grid.minor = element_blank(),
        legend.key.size  = unit(0.4, "cm"))

# Load pre-built datasets
data(FEATURE_TO_GENE)
neufer_da <- read_csv("data/neufer_da_all.csv",       show_col_types = FALSE)
mito_da   <- read_csv("data/mito_dynamics_da_all.csv", show_col_types = FALSE)
all_da    <- bind_rows(neufer_da, mito_da)       # combined for lactate genes

# ══════════════════════════════════════════════════════════════════════════════
# FOLLOW-UP 1 — DRP1 (Dnm1l) phosphoproteomics
# ══════════════════════════════════════════════════════════════════════════════
cat("=== Follow-up 1: DRP1 phosphoproteomics ===\n\n")

# Human→rat site mapping (confirmed from FEATURE_TO_GENE isoform numbering)
# NP_446107.2 canonical isoform:  S596 ≈ human S616 (CDK1/CAMKII, pro-fission)
# XP_006248777.1 isoform:         S641 ≈ human S637 (PKA, anti-fission / inhibitory)
# Additional covered sites on other isoforms annotated below
SITE_MAP <- tribble(
  ~feature_ID,                          ~site_label,        ~human_equiv,  ~known_kinase,
  "NP_446107.2_S596s",                  "S596 (rat)",       "S616",        "CDK1 / CAMKII",
  "XP_006248777.1_S641s",               "S641 (rat)",       "S637",        "PKA",
  "XP_008767063.1_S635s",               "S635 iso3 (rat)",  "~S637",       "PKA-like",
  "XP_006248784.1_S602s",               "S602 iso2 (rat)",  "~S616",       "CDK1-like",
  "XP_006248777.1_S553s",               "S553 (rat)",       "S570",        "unknown",
  "XP_006248777.1_S557s",               "S557 (rat)",       "S574",        "unknown",
  "XP_006248777.1_S563s",               "S563 (rat)",       "S580",        "unknown",
  "XP_006248777.1_S567s",               "S567 (rat)",       "S584",        "unknown",
  "XP_006248777.1_S145s",               "S145 (rat)",       "S161",        "unknown",
  "XP_008767063.1_S547s",               "S547 iso3",        "S564",        "unknown",
  "XP_008767063.1_S551s",               "S551 iso3",        "S568",        "unknown",
  "XP_008767063.1_S557s",               "S557 iso3",        "S574",        "unknown",
  "XP_008767063.1_S561s",               "S561 iso3",        "S578",        "unknown",
  "XP_008767063.1_S139s",               "S139 iso3",        "S155",        "unknown",
  "XP_006248786.1_S535s",               "S535 iso4",        "~S553",       "unknown",
  "XP_006248786.1_S585s",               "S585 iso4",        "~S602",       "unknown",
  "XP_006248786.1_S126s",               "S126 iso4",        "S143",        "unknown"
)

all_dnm1l_ids <- FEATURE_TO_GENE %>%
  filter(gene_symbol == "Dnm1l", grepl("^[NX]P_.*_S", feature_ID)) %>%
  pull(feature_ID) %>% unique()

all_items     <- data(package = "MotrpacRatTraining6moData")$results[, "Item"]
phospho_tables <- grep("^PHOSPHO_.*_DA$", all_items, value = TRUE)

cat("Loading", length(phospho_tables), "PHOSPHO tables...\n")
drp1_phospho <- map_dfr(phospho_tables, function(tbl) {
  e <- new.env()
  data(list = tbl, package = "MotrpacRatTraining6moData", envir = e)
  da <- get(tbl, envir = e)
  da %>%
    filter(feature_ID %in% all_dnm1l_ids) %>%
    select(feature_ID, assay, tissue, sex, comparison_group,
           logFC, adj_p_value, p_value,
           any_of(c("tscore","selection_fdr")))
})

drp1_phospho <- drp1_phospho %>%
  left_join(SITE_MAP, by = "feature_ID") %>%
  mutate(
    sex_label        = str_to_title(sex),
    comparison_group = factor(comparison_group, levels = TP_LEVELS),
    is_key_site      = human_equiv %in% c("S616","S637"),
    sig_dot          = if_else(adj_p_value < 0.05, "\u25cf", "")
  )

cat("DRP1 phospho rows:", nrow(drp1_phospho), "\n")
cat("Tissues covered:", paste(sort(unique(drp1_phospho$tissue)), collapse = ", "), "\n")

# Pull Dnm1l TRNSCRPT from mito data
drp1_trx <- mito_da %>%
  filter(gene_symbol == "Dnm1l", assay == "TRNSCRPT") %>%
  mutate(sex_label        = str_to_title(sex),
         comparison_group = factor(comparison_group, levels = TP_LEVELS))

# Key sites: S596 (pro-fission) and S641 (anti-fission)
key_sites <- drp1_phospho %>%
  filter(human_equiv %in% c("S616","S637"))

cat("\nKey DRP1 sites (S596/S641) across tissues:\n")
key_sites %>%
  select(site_label, human_equiv, known_kinase, tissue, sex_label,
         comparison_group, logFC, adj_p_value) %>%
  arrange(human_equiv, adj_p_value) %>%
  head(40) %>%
  as.data.frame() %>%
  print()

write_csv(drp1_phospho, "data/drp1_phospho_all_sites.csv")
write_csv(key_sites,    "data/drp1_phospho_key_sites.csv")
cat("\nSaved data/drp1_phospho_all_sites.csv\n")
cat("Saved data/drp1_phospho_key_sites.csv\n")

# --- Fig: key site trajectories + TRNSCRPT overlay per tissue ---
phospho_tissues <- intersect(unique(key_sites$tissue), unique(drp1_trx$tissue))

# Panel A: Key site logFC trajectories per tissue × sex
fig_drp1_sites <- ggplot(
  key_sites %>% filter(tissue %in% phospho_tissues),
  aes(x = comparison_group, y = logFC,
      group = interaction(feature_ID, sex_label),
      colour = sex_label, linetype = human_equiv)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey65", linewidth = 0.35) +
  geom_line(linewidth = 0.85) +
  geom_point(aes(shape = sex_label), size = 2.2) +
  geom_text(aes(label = sig_dot), nudge_y = 0.08, size = 3.5, show.legend = FALSE) +
  scale_colour_manual(values = SEX_COL, name = "Sex") +
  scale_shape_manual( values = SEX_SHAPE, name = "Sex") +
  scale_linetype_manual(
    values = c("S616" = "solid", "S637" = "dashed"),
    name   = "Human site equiv.",
    labels = c("S616" = "S596/S616 (pro-fission, CDK1/CAMKII)",
               "S637" = "S641/S637 (anti-fission, PKA)")
  ) +
  facet_wrap(~ tissue, ncol = 4) +
  labs(
    title    = "DRP1 (Dnm1l) Phosphorylation — Key Regulatory Sites",
    subtitle = "\u25cf adj_p < 0.05; solid = pro-fission (S616 equiv.), dashed = anti-fission (S637 equiv.)",
    x = "Training duration", y = "log\u2082FC (phosphosite)"
  ) +
  theme_pub +
  theme(legend.position = "bottom", legend.box = "vertical")

# Panel B: TRNSCRPT for same tissues (should be flat if PTM-only regulation)
fig_drp1_trx <- ggplot(
  drp1_trx %>% filter(tissue %in% phospho_tissues),
  aes(x = comparison_group, y = logFC,
      group = sex_label, colour = sex_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey65", linewidth = 0.35) +
  geom_line(linewidth = 0.85) +
  geom_point(aes(shape = sex_label), size = 2.2) +
  geom_text(aes(label = if_else(adj_p_value < 0.05, "\u25cf", "")),
            nudge_y = 0.05, size = 3.5, show.legend = FALSE) +
  scale_colour_manual(values = SEX_COL, name = "Sex") +
  scale_shape_manual( values = SEX_SHAPE, name = "Sex") +
  facet_wrap(~ tissue, ncol = 4) +
  labs(
    title    = "DRP1 (Dnm1l) Transcript — Reference for PTM vs. Abundance",
    subtitle = "Compare to phosphosite panel above; flat transcript = pure PTM regulation",
    x = "Training duration", y = "log\u2082FC (TRNSCRPT)"
  ) +
  theme_pub + theme(legend.position = "bottom")

fig_drp1 <- fig_drp1_sites / fig_drp1_trx +
  plot_annotation(title = "DRP1 Post-translational Regulation of Fission",
                  tag_levels = "A")

ggsave("figures/fig_followup_drp1_phospho.pdf",
       fig_drp1, width = 14, height = 12, device = cairo_pdf)
ggsave("figures/fig_followup_drp1_phospho.png",
       fig_drp1, width = 14, height = 12, dpi = 180)
cat("Saved figures/fig_followup_drp1_phospho\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# FOLLOW-UP 2 — Lactate axis tissue coordination (producer vs. consumer map)
# ══════════════════════════════════════════════════════════════════════════════
cat("=== Follow-up 2: Lactate axis tissue coordination ===\n\n")

# Gene roles:
#   EXPORT (glycolytic / lactate producing): Ldha + Slc16a3 (MCT4)
#   IMPORT (oxidative / lactate consuming):  Ldhb + Slc16a1 (MCT1)
LACTATE_GENES <- c("Slc16a1","Slc16a3","Ldha","Ldhb")
GENE_ROLE     <- c(Slc16a1 = "Import (MCT1)", Slc16a3 = "Export (MCT4)",
                   Ldha    = "Export (LDHA)", Ldhb    = "Import (LDHB)")
GENE_ORDER    <- c("Slc16a3","Ldha","Slc16a1","Ldhb")   # export → import
ROLE_COL      <- c("Export (MCT4)" = "#D6604D", "Export (LDHA)" = "#F4A582",
                   "Import (MCT1)" = "#2166AC", "Import (LDHB)" = "#92C5DE")

lac_all <- all_da %>%
  filter(gene_symbol %in% LACTATE_GENES, assay == "TRNSCRPT") %>%
  mutate(
    sex_label        = str_to_title(sex),
    comparison_group = factor(comparison_group, levels = TP_LEVELS),
    tissue           = factor(tissue, levels = TISSUE_ORDER),
    gene_symbol      = factor(gene_symbol, levels = GENE_ORDER),
    role             = GENE_ROLE[as.character(gene_symbol)],
    sig_dot          = if_else(adj_p_value < 0.05, "\u25cf", "")
  ) %>%
  distinct(gene_symbol, tissue, sex_label, comparison_group, .keep_all = TRUE)

# ETS index at 8w per tissue × sex (for overlay)
ets_8w <- neufer_da %>%
  filter(assay == "TRNSCRPT", comparison_group == "8w",
         framework_category %in% c("ets_CI","ets_CIII","ets_CIV")) %>%
  group_by(tissue, sex) %>%
  summarise(ets_mean = mean(logFC, na.rm = TRUE),
            ets_n_sig = sum(adj_p_value < 0.05, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(sex_label = str_to_title(sex),
         tissue    = factor(tissue, levels = TISSUE_ORDER))

lac_8w <- lac_all %>% filter(comparison_group == "8w")

# Panel A: tissue × gene heatmap at 8w (both sexes)
lim_lac <- ceiling(max(abs(lac_8w$logFC), na.rm = TRUE) * 10) / 10

fig_lac_heat <- ggplot(lac_8w,
                       aes(x = gene_symbol, y = tissue, fill = logFC)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_dot), size = 2.8, colour = "grey15", vjust = 0.75) +
  scale_fill_distiller(palette = DIVPAL, limits = c(-lim_lac, lim_lac),
                       direction = -1, na.value = "grey88",
                       name = "log\u2082FC\n(TRNSCRPT, 8w)") +
  scale_y_discrete(limits = rev(TISSUE_ORDER[TISSUE_ORDER %in% unique(lac_8w$tissue)])) +
  scale_x_discrete(labels = c(Slc16a3 = "Slc16a3\n(MCT4 export)",
                               Ldha    = "Ldha\n(LDHA export)",
                               Slc16a1 = "Slc16a1\n(MCT1 import)",
                               Ldhb    = "Ldhb\n(LDHB import)")) +
  # Separator between export and import genes
  geom_vline(xintercept = 2.5, colour = "grey30", linewidth = 0.8) +
  annotate("text", x = 1.5, y = 20.3, label = "\u2190 Export / Glycolytic",
           size = 2.8, fontface = "bold", hjust = 0.5) +
  annotate("text", x = 3.5, y = 20.3, label = "Import / Oxidative \u2192",
           size = 2.8, fontface = "bold", hjust = 0.5) +
  facet_wrap(~ sex_label, ncol = 2) +
  coord_cartesian(clip = "off") +
  labs(title    = "Lactate Axis — Tissue × Gene at 8w (TRNSCRPT)",
       subtitle = "\u25cf adj_p < 0.05 | Left = lactate producers, Right = lactate consumers",
       x = NULL, y = NULL) +
  theme_pub +
  theme(axis.text.x  = element_text(size = 8.5),
        axis.text.y  = element_text(size = 8),
        panel.border = element_blank(),
        plot.margin  = margin(t = 22, r = 8, b = 5, l = 5))

# Panel B: Lactate import/export index vs ETS index across tissues
# Export index = mean(Ldha, Slc16a3) at 8w; Import index = mean(Ldhb, Slc16a1)
lac_index <- lac_8w %>%
  mutate(axis = if_else(gene_symbol %in% c("Ldha","Slc16a3"), "export", "import")) %>%
  group_by(tissue, sex_label, axis) %>%
  summarise(mean_lfc = mean(logFC, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = axis, values_from = mean_lfc) %>%
  left_join(ets_8w, by = c("tissue","sex_label")) %>%
  mutate(net_import = import - export)   # positive = net import orientation

fig_lac_ets <- ggplot(lac_index,
                      aes(x = ets_mean, y = net_import,
                          colour = sex_label, shape = sex_label, label = tissue)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.12,
              aes(x = ets_mean, y = net_import, colour = sex_label, fill = sex_label),
              inherit.aes = FALSE) +
  geom_point(size = 3) +
  geom_text_repel(size = 2.6, max.overlaps = 18, seed = 42,
                  segment.size = 0.25, show.legend = FALSE) +
  scale_colour_manual(values = SEX_COL, name = "Sex") +
  scale_shape_manual( values = SEX_SHAPE, name = "Sex") +
  scale_fill_manual(  values = SEX_COL, name = "Sex") +
  annotate("text", x = -Inf, y = Inf, label = "Net importer\n(Ldhb + MCT1 > Ldha + MCT4)",
           hjust = -0.05, vjust = 1.3, size = 2.8, colour = "grey30") +
  annotate("text", x = -Inf, y = -Inf, label = "Net exporter\n(Ldha + MCT4 > Ldhb + MCT1)",
           hjust = -0.05, vjust = -0.3, size = 2.8, colour = "grey30") +
  facet_wrap(~ sex_label) +
  labs(
    title    = "Lactate Import/Export Balance vs. ETS Expansion Across Tissues (8w)",
    subtitle = "Net import index = mean(Slc16a1 + Ldhb) \u2212 mean(Slc16a3 + Ldha) log\u2082FC",
    x        = "ETS mean log\u2082FC (CI+III+IV, TRNSCRPT, 8w)",
    y        = "Net lactate import orientation (log\u2082FC)"
  ) +
  theme_pub + theme(legend.position = "none")

# Panel C: Temporal trajectories for significant genes
lac_sig_genes <- lac_all %>%
  group_by(gene_symbol, tissue, sex_label) %>%
  filter(any(adj_p_value < 0.05)) %>%
  ungroup()

if (nrow(lac_sig_genes) > 0) {
  fig_lac_traj <- ggplot(lac_sig_genes,
                         aes(x = comparison_group, y = logFC,
                             group = interaction(tissue, sex_label),
                             colour = tissue, linetype = sex_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey65", linewidth = 0.35) +
    geom_line(linewidth = 0.9) +
    geom_point(aes(shape = sex_label), size = 2.2) +
    geom_text(aes(label = sig_dot), nudge_x = 0.2, size = 3, show.legend = FALSE) +
    scale_colour_brewer(palette = "Dark2", name = "Tissue") +
    scale_linetype_manual(values = c(Male = "solid", Female = "dashed"), name = "Sex") +
    scale_shape_manual(values = SEX_SHAPE, name = "Sex") +
    facet_wrap(~ gene_symbol, nrow = 1,
               labeller = labeller(gene_symbol = GENE_ROLE)) +
    labs(
      title    = "Lactate Axis — Temporal Trajectories (significant tissue × gene pairs)",
      subtitle = "\u25cf adj_p < 0.05",
      x = "Training duration", y = "log\u2082FC"
    ) +
    theme_pub + theme(legend.position = "bottom", legend.box = "horizontal")
} else {
  fig_lac_traj <- ggplot() + theme_void() +
    labs(title = "No significant temporal trajectories in lactate genes")
}

fig_lac <- (fig_lac_heat / (fig_lac_ets | fig_lac_traj)) +
  plot_annotation(title = "Lactate Axis Tissue Coordination",
                  tag_levels = "A") +
  plot_layout(heights = c(1.4, 1))

ggsave("figures/fig_followup_lactate_coordination.pdf",
       fig_lac, width = 14, height = 14, device = cairo_pdf)
ggsave("figures/fig_followup_lactate_coordination.png",
       fig_lac, width = 14, height = 14, dpi = 180)
cat("Saved figures/fig_followup_lactate_coordination\n")

write_csv(lac_all,   "data/lactate_axis_da_all.csv")
write_csv(lac_index, "data/lactate_import_export_index.csv")
cat("Saved data/lactate_axis_da_all.csv\n")
cat("Saved data/lactate_import_export_index.csv\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# FOLLOW-UP 3 — Succinate-SDH connection in lung
# ══════════════════════════════════════════════════════════════════════════════
cat("=== Follow-up 3: Succinate-SDH in lung ===\n\n")

# SDH subunits from neufer analysis (LUNG, all timepoints, TRNSCRPT + PROT)
sdh_lung <- neufer_da %>%
  filter(gene_symbol %in% c("Sdha","Sdhb","Sdhc","Sdhd"),
         tissue == "LUNG") %>%
  mutate(
    sex_label        = str_to_title(sex),
    comparison_group = factor(comparison_group, levels = TP_LEVELS),
    sig_dot          = if_else(adj_p_value < 0.05, "\u25cf", ""),
    assay_label      = case_when(
      assay == "TRNSCRPT" ~ "Transcript",
      assay == "PROT"     ~ "Protein",
      TRUE                ~ assay
    )
  )

cat("SDH subunits in LUNG (TRNSCRPT + PROT):\n")
print(sdh_lung %>%
      select(gene_symbol, assay, sex_label, comparison_group, logFC, adj_p_value) %>%
      arrange(gene_symbol, assay, sex_label, comparison_group), n = 60)

# Succinate and fumarate metabolites in LUNG
data(METAB_LUNG_DA)
lung_metab <- METAB_LUNG_DA %>%
  filter(feature_ID %in% c("succinate","Succinate","fumarate","Fumarate",
                            "malate","Malate","citrate","Citrate",
                            "NAD+","NADH","NADP+","NADPH","FAD","pyruvate","lactate")) %>%
  mutate(
    canonical = case_when(
      grepl("^[Ss]uccinate$",  feature_ID) ~ "Succinate",
      grepl("^[Ff]umarate$",   feature_ID) ~ "Fumarate",
      grepl("^[Mm]alate$",     feature_ID) ~ "Malate",
      grepl("^[Cc]itrate$",    feature_ID) ~ "Citrate",
      grepl("^[Ll]actate$",    feature_ID) ~ "Lactate",
      grepl("^[Pp]yruvate$",   feature_ID) ~ "Pyruvate",
      TRUE ~ feature_ID
    ),
    sex_label        = str_to_title(sex),
    comparison_group = factor(comparison_group, levels = TP_LEVELS),
    sig_dot          = if_else(adj_p_value < 0.05, "\u25cf", "")
  ) %>%
  group_by(canonical, sex_label, comparison_group) %>%
  summarise(logFC = mean(logFC, na.rm = TRUE),
            adj_p_value = min(adj_p_value, na.rm = TRUE),
            p_value     = min(p_value, na.rm = TRUE),
            sig_dot = if_else(min(adj_p_value, na.rm = TRUE) < 0.05, "\u25cf", ""),
            .groups = "drop")

cat("\nLUNG metabolites (all timepoints):\n")
print(lung_metab %>% arrange(canonical, sex_label, comparison_group), n = 60)

# SUCNR1 status note
cat("\nNote: Sucnr1 (GPR91) is present in FEATURE_TO_GENE (ENSRNOG00000014039) but\n")
cat("has zero rows in TRNSCRPT_LUNG_DA — gene not detected/expressed in lung tissue.\n\n")

# SDH mean per timepoint (mean of Sdha/b/c/d, separate by assay and sex)
sdh_mean <- sdh_lung %>%
  group_by(assay_label, sex_label, comparison_group) %>%
  summarise(
    mean_lfc = mean(logFC, na.rm = TRUE),
    se_lfc   = sd(logFC, na.rm = TRUE) / sqrt(n_distinct(gene_symbol)),
    n_sig    = sum(adj_p_value < 0.05, na.rm = TRUE),
    .groups  = "drop"
  )

# Succinate metabolite only for the primary plot
succinate_lung <- lung_metab %>%
  filter(canonical == "Succinate")

cat("Succinate in LUNG:\n")
print(succinate_lung)

# --- Fig: dual-panel — SDH genes (top) + succinate metabolite (bottom) ---

# Panel A: SDH transcript + protein trajectories (all 4 subunits)
fig_sdh_genes <- ggplot(sdh_lung,
                        aes(x = comparison_group, y = logFC,
                            group = interaction(gene_symbol, sex_label),
                            colour = gene_symbol, linetype = sex_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey65", linewidth = 0.35) +
  geom_ribbon(
    data = sdh_mean,
    aes(x = comparison_group, ymin = mean_lfc - se_lfc, ymax = mean_lfc + se_lfc,
        group = interaction(assay_label, sex_label)),
    fill = "grey70", alpha = 0.2, colour = NA, inherit.aes = FALSE
  ) +
  geom_line(aes(group = interaction(gene_symbol, assay_label, sex_label)),
            linewidth = 0.8) +
  geom_point(aes(shape = sex_label), size = 2) +
  geom_text(aes(label = sig_dot), nudge_y = 0.03, size = 3.5, show.legend = FALSE) +
  scale_colour_brewer(palette = "Set1", name = "SDH subunit") +
  scale_linetype_manual(values = c(Male = "solid", Female = "dashed"), name = "Sex") +
  scale_shape_manual(values = SEX_SHAPE, name = "Sex") +
  facet_wrap(~ assay_label, ncol = 2) +
  labs(
    title    = "SDH Subunits in LUNG — Transcript and Protein Trajectories",
    subtitle = "Shaded band = mean \u00b1 SE across Sdha/b/c/d; \u25cf adj_p < 0.05",
    x = "Training duration", y = "log\u2082FC vs. sedentary"
  ) +
  theme_pub + theme(legend.position = "right")

# Panel B: Succinate metabolite trajectory (LUNG)
fig_succinate <- ggplot(succinate_lung,
                        aes(x = comparison_group, y = logFC,
                            group = sex_label, colour = sex_label, shape = sex_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey65", linewidth = 0.35) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  geom_text(aes(label = sig_dot), nudge_y = 0.05, size = 4, show.legend = FALSE) +
  scale_colour_manual(values = SEX_COL, name = "Sex") +
  scale_shape_manual( values = SEX_SHAPE, name = "Sex") +
  labs(
    title    = "Succinate Metabolite in LUNG",
    subtitle = "\u25cf adj_p < 0.05 — note: SDH subunits are flat above",
    x = "Training duration", y = "log\u2082FC vs. sedentary"
  ) +
  theme_pub + theme(legend.position = "right")

# Panel C: All TCA metabolites in LUNG for context
tca_metabs_lung <- lung_metab %>%
  filter(canonical %in% c("Succinate","Fumarate","Malate","Citrate","Lactate","Pyruvate"))

fig_lung_tca <- ggplot(tca_metabs_lung,
                       aes(x = comparison_group, y = logFC,
                           group = interaction(canonical, sex_label),
                           colour = canonical, linetype = sex_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey65", linewidth = 0.35) +
  geom_line(linewidth = 0.85) +
  geom_point(aes(shape = sex_label), size = 2) +
  geom_text(aes(label = sig_dot), nudge_x = 0.2, size = 3.5, show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2", name = "Metabolite") +
  scale_linetype_manual(values = c(Male = "solid", Female = "dashed"), name = "Sex") +
  scale_shape_manual(values = SEX_SHAPE, name = "Sex") +
  labs(
    title    = "TCA / Glycolytic Metabolites in LUNG",
    subtitle = "\u25cf adj_p < 0.05",
    x = "Training duration", y = "log\u2082FC vs. sedentary"
  ) +
  theme_pub + theme(legend.position = "right")

fig_sdh <- (fig_sdh_genes / (fig_succinate | fig_lung_tca)) +
  plot_annotation(
    title    = "Succinate-SDH Discordance in Lung: Metabolite Down, Enzyme Flat",
    subtitle = "Succinate significantly decreases across timepoints without corresponding SDH regulation",
    tag_levels = "A"
  ) +
  plot_layout(heights = c(1, 0.85))

ggsave("figures/fig_followup_lung_succinate_sdh.pdf",
       fig_sdh, width = 13, height = 12, device = cairo_pdf)
ggsave("figures/fig_followup_lung_succinate_sdh.png",
       fig_sdh, width = 13, height = 12, dpi = 180)
cat("Saved figures/fig_followup_lung_succinate_sdh\n")

write_csv(sdh_lung,     "data/sdh_lung_all.csv")
write_csv(lung_metab,   "data/lung_tca_metabolites.csv")
cat("Saved data/sdh_lung_all.csv\n")
cat("Saved data/lung_tca_metabolites.csv\n")

# ── Final summary ─────────────────────────────────────────────────────────────
cat("\n=== SUMMARY OF FOLLOW-UP FINDINGS ===\n\n")

cat("1. DRP1 phosphoproteomics:\n")
sig_drp1 <- drp1_phospho %>% filter(adj_p_value < 0.05)
cat("   Significant phosphosites:", nrow(sig_drp1), "\n")
if (nrow(sig_drp1) > 0)
  print(sig_drp1 %>% select(site_label, human_equiv, tissue, sex_label,
                             comparison_group, logFC, adj_p_value) %>%
        arrange(adj_p_value), n = 20)

cat("\n2. Lactate axis significant genes:\n")
print(lac_all %>% filter(adj_p_value < 0.05) %>%
      select(gene_symbol, tissue, sex_label, comparison_group, logFC, adj_p_value) %>%
      arrange(adj_p_value), n = 20)

cat("\n3. Lung succinate significant timepoints:\n")
print(succinate_lung %>% filter(adj_p_value < 0.05) %>%
      select(canonical, sex_label, comparison_group, logFC, adj_p_value))

cat("\n   SDH subunit significant hits in LUNG: ",
    sum(sdh_lung$adj_p_value < 0.05), "(should be 0 — enzyme flat)\n")

cat("\nAll files saved. Done.\n")
