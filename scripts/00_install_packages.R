# =============================================================================
# 00_install_packages.R
# Purpose : Install all R packages required to reproduce this analysis.
# Inputs  : none
# Outputs : none (packages installed to local R library)
# Run once before any other script.
# =============================================================================

# ── Bioconductor packages (MoTrPAC data) ─────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

BiocManager::install(
  c("MotrpacRatTraining6moData",   # differential analysis results + normalized data
    "MotrpacRatTraining6mo",        # analysis helper functions
    "ComplexHeatmap",               # gene-by-gene correlation heatmaps
    "circlize"),                    # color scale helpers for ComplexHeatmap
  ask = FALSE, update = FALSE
)

# ── CRAN packages ─────────────────────────────────────────────────────────────
cran_pkgs <- c(
  "tidyverse",   # data wrangling + ggplot2
  "patchwork",   # multi-panel figure layout
  "ggrepel",     # non-overlapping text labels
  "broom",       # tidy model summaries
  "readxl"       # reading human supplement XLSX files (07_human_data.R)
)

install.packages(
  cran_pkgs[!cran_pkgs %in% rownames(installed.packages())],
  repos = "https://cloud.r-project.org"
)

# ── Python packages (sensitivity table PDF) ───────────────────────────────────
# Used by supp_sensitivity_table_pdf.py — run separately if needed:
#   pip3 install reportlab Pillow

message("\nAll packages installed. Run scripts in order: 01 → 02 → 03 → 04 → 05 → 06 → 07")
