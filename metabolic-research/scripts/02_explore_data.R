# 02_explore_data.R — inspect raw data structure before analysis

source(file.path("scripts", "config.R"))
library(tidyverse)

raw_files <- list.files(DIR_DATA_RAW, full.names = TRUE, recursive = TRUE)
raw_files <- raw_files[!grepl("\\.gitkeep$", raw_files)]

if (length(raw_files) == 0) {
  stop("No raw data found. Run 01_download_data.R first.")
}

cat("Raw files:\n"); print(raw_files)

# Inspect dimensions, column types, missingness, summary stats here.
