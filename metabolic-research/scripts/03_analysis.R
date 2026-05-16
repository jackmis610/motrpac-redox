# 03_analysis.R — core reanalysis
#
# Read from data/raw/ (or data/processed/), run the analysis, and write
# tables to results/ and plots to figures/.

source(file.path("scripts", "config.R"))
library(tidyverse)

# --- Load --------------------------------------------------------------------

# --- Analyze -----------------------------------------------------------------

# --- Save --------------------------------------------------------------------
# write_csv(result_table, file.path(DIR_RESULTS, "main_result.csv"))
# ggsave(file.path(DIR_FIGURES, "main_figure.pdf"), plot = p,
#        device = PDF_DEVICE, width = 7, height = 5)
