## ============================================================
## Core package installer for Code Ocean
## ============================================================
## Run this file only if the Code Ocean R environment is missing
## core CRAN packages required by the default reviewer-safe workflow.
##
## Recommended order in a fresh Code Ocean capsule:
##   1. bash codeocean_system_deps.sh        # only if system libraries are needed
##   2. Rscript 00_install_core_packages.R   # install core CRAN packages
##   3. Rscript packages.R                   # check installed packages
##   4. Rscript run_codeocean.R              # run reviewer-safe workflow
##
## This script intentionally does not install heavy optional packages:
##   monocle3, infercnv, copykat, SCENIC, AUCell, RcisTarget,
##   GENIE3, GSVA, ComplexHeatmap, Nebulosa
## These packages are skipped by run_codeocean.R unless available.
## ============================================================

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

dir.create("output", showWarnings = FALSE, recursive = TRUE)

core_cran_packages <- c(
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

install_one <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    message(pkg, " is already installed.")
    return(TRUE)
  }

  message("\nInstalling core CRAN package: ", pkg)

  ok <- tryCatch(
    {
      install.packages(
        pkg,
        repos = "https://cloud.r-project.org",
        dependencies = c("Depends", "Imports", "LinkingTo")
      )
      requireNamespace(pkg, quietly = TRUE)
    },
    error = function(e) {
      message("Failed to install ", pkg, ": ", conditionMessage(e))
      FALSE
    },
    warning = function(w) {
      message("Warning while installing ", pkg, ": ", conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  isTRUE(ok)
}

status <- data.frame(
  package = core_cran_packages,
  installed_before = vapply(core_cran_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)),
  installed_after = NA,
  version = NA_character_,
  stringsAsFactors = FALSE
)

for (pkg in core_cran_packages) {
  install_one(pkg)
}

status$installed_after <- vapply(core_cran_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))

for (pkg in status$package[status$installed_after]) {
  status$version[status$package == pkg] <- as.character(utils::packageVersion(pkg))
}

utils::write.csv(
  status,
  file = file.path("output", "core_package_install_check.csv"),
  row.names = FALSE
)

missing_final <- status$package[!status$installed_after]

if (length(missing_final) > 0) {
  message("\nThe following core packages are still missing:")
  message(paste(missing_final, collapse = ", "))
  message("Check the package installation log and install missing packages through the Code Ocean environment if needed.")
} else {
  message("\nAll core CRAN packages are installed.")
}

message("\nCore package installation check written to output/core_package_install_check.csv")
