## ============================================================
## Package checklist for YB_1_scRNA
## ============================================================
## This file checks whether the main R packages used by the scripts are
## installed. It does not force installation of Bioconductor/GitHub-only
## packages, because those packages often require version-specific setup.
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
  "Nebulosa",
  "pheatmap",
  "ggpubr",
  "cowplot"
)

special_packages <- c(
  "monocle3",
  "infercnv",
  "copykat",
  "SCENIC",
  "AUCell",
  "RcisTarget",
  "GENIE3",
  "GSVA",
  "ComplexHeatmap"
)

all_packages <- unique(c(cran_packages, special_packages))
installed <- rownames(installed.packages())
missing <- setdiff(all_packages, installed)

if (length(missing) == 0) {
  message("All listed packages are installed.")
} else {
  message("Missing packages:")
  message(paste(missing, collapse = ", "))
  message("\nInstall missing CRAN packages with install.packages().")
  message("Install Bioconductor/GitHub packages according to their package documentation.")
}

## Optional CRAN installation helper. Uncomment if needed.
# missing_cran <- setdiff(cran_packages, rownames(installed.packages()))
# if (length(missing_cran) > 0) {
#   install.packages(missing_cran)
# }
