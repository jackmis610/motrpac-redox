# 00_install_packages.R — install dependencies (run once)

cran_packages <- c(
  "rprojroot",   # robust project-root detection (used by config.R)
  "tidyverse",   # data wrangling + ggplot2
  "data.table",  # fast IO for large tables
  "here"         # path helper
)

to_install <- cran_packages[!cran_packages %in% rownames(installed.packages())]
if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

# Bioconductor packages — uncomment and add as your reanalysis needs them.
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager", repos = "https://cloud.r-project.org")
# }
# BiocManager::install(c("GEOquery", "limma", "DESeq2"))

cat("Package setup complete.\n")
