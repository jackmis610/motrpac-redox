# 01_download_data.R — fetch public datasets into data/raw/
#
# Keep every download scripted and reproducible. Raw files in data/raw/ are
# gitignored, so anyone can rebuild them by re-running this script.

source(file.path("scripts", "config.R"))

# --- Example: download a file from a public URL ---------------------------
# dataset_url  <- "https://example.org/path/to/dataset.csv"
# dest         <- file.path(DIR_DATA_RAW, "dataset.csv")
# if (!file.exists(dest)) {
#   download.file(dataset_url, dest, mode = "wb")
# }

# --- Example: pull a GEO series with GEOquery ------------------------------
# library(GEOquery)
# gse <- getGEO("GSEXXXXXX", destdir = DIR_DATA_RAW)

cat("Add dataset download steps above, then re-run this script.\n")
