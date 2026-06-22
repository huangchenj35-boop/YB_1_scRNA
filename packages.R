## ============================================================
## Package checklist / optional installer for YB_1_scRNA
## ============================================================
## This file separates CRAN, Bioconductor, and GitHub packages.
## Nebulosa is not a CRAN package, so it should not be installed with
## install.packages().
##
## Default behavior:
##   - check missing packages
##   - do not install automatically
##
## To install available CRAN packages:
##   install_cran <- TRUE
##
## To install Bioconductor packages:
##   install_bioc <- TRUE
##
## GitHub packages usually depend on the R/Bioconductor version and are
## listed separately for manual installation if needed.
## ============================================================

install_cran <- FALSE
install_bioc <- FALSE

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

github_or_manual_packages <- c(
  "monocle3",
  "copykat",
  "SCENIC"
)

all_packages <- unique(c(cran_packages, bioc_packages, github_or_manual_packages))
installed <- rownames(installed.packages())
missing <- setdiff(all_packages, installed)

if (length(missing) == 0) {
  message("All listed packages are installed.")
} else {
  message("Missing packages:")
  message(paste(missing, collapse = ", "))
}

missing_cran <- setdiff(cran_packages, installed)
missing_bioc <- setdiff(bioc_packages, installed)
missing_manual <- setdiff(github_or_manual_packages, installed)

if (length(missing_cran) > 0) {
  message("\nMissing CRAN packages:")
  message(paste(missing_cran, collapse = ", "))
  if (isTRUE(install_cran)) {
    install.packages(missing_cran, repos = "https://cloud.r-project.org")
  } else {
    message("Set install_cran <- TRUE to install these CRAN packages.")
  }
}

if (length(missing_bioc) > 0) {
  message("\nMissing Bioconductor packages:")
  message(paste(missing_bioc, collapse = ", "))
  if (isTRUE(install_bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  } else {
    message("Set install_bioc <- TRUE to install these Bioconductor packages.")
  }
}

if (length(missing_manual) > 0) {
  message("\nPackages requiring GitHub/manual installation:")
  message(paste(missing_manual, collapse = ", "))
  message("Install these according to the package-specific instructions and the R version used in Code Ocean.")
}
