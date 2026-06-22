## ============================================================
## Optional package installer for Code Ocean
## ============================================================
## This script installs common CRAN and Bioconductor packages used by
## the workflow. GitHub/manual packages are listed but not installed here.
## ============================================================

cran_packages <- c(
  "Seurat",
  "dplyr",
  "tidyr",
  "readr",
  "stringr",
  "ggplot2",
  "patchwork",
  "Matrix",
  "scales",
  "colorspace",
  "aplot",
  "pheatmap",
  "ggpubr",
  "cowplot"
)

bioc_packages <- c(
  "Nebulosa",
  "infercnv",
  "AUCell",
  "RcisTarget",
  "GENIE3",
  "GSVA",
  "ComplexHeatmap"
)

installed <- rownames(installed.packages())

missing_cran <- setdiff(cran_packages, installed)
if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

installed <- rownames(installed.packages())
missing_bioc <- setdiff(bioc_packages, installed)
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

message("Package installation check finished.")
message("Manual/GitHub packages to check separately: monocle3, copykat, SCENIC")
