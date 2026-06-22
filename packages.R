## ============================================================
## Package checklist for YB_1_scRNA
## ============================================================
## This file checks package availability only.
## It does not install packages at runtime.
##
## Recommended Code Ocean order:
##   1. bash codeocean_system_deps.sh        # optional system libraries
##   2. Rscript 00_install_core_packages.R   # install core CRAN packages
##   3. Rscript packages.R                   # check all listed packages
##   4. Rscript run_codeocean.R              # reviewer-safe workflow
##
## Heavy steps such as inferCNV, CopyKAT, SCENIC, Monocle3, and
## oncoPredict are handled as optional steps by run_codeocean.R if the
## corresponding packages or input objects are not available.
## ============================================================

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

core_packages <- c(
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

bioconductor_packages <- c(
  "Nebulosa",
  "AUCell",
  "RcisTarget",
  "GENIE3",
  "GSVA",
  "ComplexHeatmap"
)

heavy_or_manual_packages <- c(
  "monocle3",
  "infercnv",
  "copykat",
  "SCENIC",
  "oncoPredict"
)

all_packages <- unique(c(core_packages, bioconductor_packages, heavy_or_manual_packages))

package_status <- data.frame(
  package = all_packages,
  group = c(
    rep("core", length(core_packages)),
    rep("bioconductor_or_optional", length(bioconductor_packages)),
    rep("heavy_or_manual", length(heavy_or_manual_packages))
  ),
  installed = vapply(all_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)),
  version = NA_character_,
  stringsAsFactors = FALSE
)

for (pkg in package_status$package[package_status$installed]) {
  package_status$version[package_status$package == pkg] <- as.character(utils::packageVersion(pkg))
}

print(package_status, row.names = FALSE)

dir.create("output", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(
  package_status,
  file = file.path("output", "package_check.csv"),
  row.names = FALSE
)

missing_core <- package_status$package[package_status$group == "core" & !package_status$installed]
missing_optional <- package_status$package[package_status$group != "core" & !package_status$installed]

if (length(missing_core) > 0) {
  message("\nMissing core packages:")
  message(paste(missing_core, collapse = ", "))
  message("Run Rscript 00_install_core_packages.R, then run this check again.")
}

if (length(missing_optional) > 0) {
  message("\nMissing optional/heavy packages:")
  message(paste(missing_optional, collapse = ", "))
  message("run_codeocean.R will skip steps that require these packages unless they are installed.")
}

if (length(missing_core) == 0 && length(missing_optional) == 0) {
  message("\nAll listed packages are installed.")
} else {
  message("\nPackage check finished with missing packages. See output/package_check.csv.")
}
