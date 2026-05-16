# config.R — shared paths and settings
# Sourced at the top of every numbered script.

# Resolve project root regardless of where R is launched from.
PROJECT_ROOT <- tryCatch(
  rprojroot::find_root(rprojroot::has_file(".gitignore")),
  error = function(e) getwd()
)

DIR_DATA_RAW       <- file.path(PROJECT_ROOT, "data", "raw")
DIR_DATA_PROCESSED <- file.path(PROJECT_ROOT, "data", "processed")
DIR_FIGURES        <- file.path(PROJECT_ROOT, "figures")
DIR_RESULTS        <- file.path(PROJECT_ROOT, "results")

for (d in c(DIR_DATA_RAW, DIR_DATA_PROCESSED, DIR_FIGURES, DIR_RESULTS)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Reproducibility
set.seed(1)

# PDF output that handles Unicode characters cleanly.
PDF_DEVICE <- grDevices::cairo_pdf
