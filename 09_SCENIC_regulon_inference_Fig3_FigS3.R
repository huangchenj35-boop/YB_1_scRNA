## ============================================================
## SCENIC Regulon Inference
## Figures: data source for Fig 3A, Fig 3B, Fig S3A, Fig S3B, Fig S3D
## Dataset: GSE138709 (5 tumor + 3 adjacent iCCA samples)
## ============================================================
## Prerequisites:
##   - scRNA1: Seurat object from 01_Fig1A_FigS1B_cell_annotation.R
##   - SCENIC motif database feather files in scenic_db_dir / SCENIC_DB
## Output:
##   - int/3.4_regulonAUC.Rds
##   - int/4.1_binaryRegulonActivity.Rds
##   - output/SCENIC_regulon_target_table.tsv, when available
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SCENIC)
  library(AUCell)
  library(RcisTarget)
  library(GENIE3)
})

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
set.seed(1234)

if (!exists("scRNA1") || !inherits(scRNA1, "Seurat")) {
  candidates <- c("output/scRNA1_annotated.rds", "output/scRNA1_preprocessed.rds")
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) {
    stop("scRNA1 was not found in memory and no expected RDS file exists in output/.")
  }
  scRNA1 <- readRDS(hit)
}

dir.create("int", showWarnings = FALSE, recursive = TRUE)
dir.create("output", showWarnings = FALSE, recursive = TRUE)

if (!exists("scenic_db_dir")) scenic_db_dir <- "SCENIC_DB"
if (!dir.exists(scenic_db_dir)) {
  stop("SCENIC database directory was not found: ", scenic_db_dir)
}

if (!exists("scenic_dbs")) {
  scenic_dbs <- c(
    "hg19-500bp-upstream-7species.mc9nr.feather",
    "hg19-tss-centered-10kb-7species.mc9nr.feather"
  )
  if (!all(file.exists(file.path(scenic_db_dir, scenic_dbs)))) {
    scenic_dbs <- list.files(scenic_db_dir, pattern = "\\.feather$", full.names = FALSE)
  }
}
if (!length(scenic_dbs)) stop("No SCENIC feather database files found in ", scenic_db_dir)
names(scenic_dbs) <- make.names(sub("\\.feather$", "", scenic_dbs))

DefaultAssay(scRNA1) <- if ("RNA" %in% names(scRNA1@assays)) "RNA" else DefaultAssay(scRNA1)
expr_mat <- as.matrix(GetAssayData(scRNA1, assay = DefaultAssay(scRNA1), slot = "counts"))
expr_mat <- expr_mat[rowSums(expr_mat) > 0, , drop = FALSE]

scenic_options_file <- file.path("int", "scenicOptions.Rds")
if (file.exists(scenic_options_file)) {
  scenicOptions <- readRDS(scenic_options_file)
} else {
  scenicOptions <- initializeScenic(
    org = "hgnc",
    dbDir = scenic_db_dir,
    dbs = scenic_dbs,
    datasetTitle = "GSE138709_iCCA_SCENIC",
    nCores = if (exists("nCores")) nCores else 4
  )
  saveRDS(scenicOptions, scenic_options_file)
}

genes_kept_file <- file.path("int", "exprMat_filtered_for_SCENIC.rds")
if (file.exists(genes_kept_file)) {
  expr_mat_filtered <- readRDS(genes_kept_file)
} else {
  genes_kept <- geneFiltering(expr_mat, scenicOptions)
  expr_mat_filtered <- expr_mat[genes_kept, , drop = FALSE]
  saveRDS(expr_mat_filtered, genes_kept_file)
}

if (!file.exists(file.path("int", "1.1_genie3ll.Rds"))) {
  runCorrelation(expr_mat_filtered, scenicOptions)
  runGenie3(log2(expr_mat_filtered + 1), scenicOptions)
}

if (!file.exists(file.path("int", "1.4_GENIE3_linkList.Rds"))) {
  scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
  saveRDS(scenicOptions, scenic_options_file)
}

if (!file.exists(file.path("int", "2.5_regulonTargetsInfo.Rds"))) {
  scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
  saveRDS(scenicOptions, scenic_options_file)
}

if (!file.exists(file.path("int", "3.4_regulonAUC.Rds"))) {
  scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat = expr_mat)
  saveRDS(scenicOptions, scenic_options_file)
}

if (!file.exists(file.path("int", "4.1_binaryRegulonActivity.Rds"))) {
  scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
  saveRDS(scenicOptions, scenic_options_file)
}

target_rds <- file.path("int", "2.5_regulonTargetsInfo.Rds")
if (file.exists(target_rds)) {
  targets <- readRDS(target_rds)
  if (is.data.frame(targets)) {
    data.table::fwrite(targets, "output/SCENIC_regulon_target_table.tsv", sep = "\t")
  }
}

saveRDS(scenicOptions, "output/scenicOptions.rds")
message("Done: SCENIC regulon inference completed or reused from int/.")
