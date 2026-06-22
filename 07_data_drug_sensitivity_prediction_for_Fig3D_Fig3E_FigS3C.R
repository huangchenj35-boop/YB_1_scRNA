## ============================================================
## Drug Sensitivity Prediction (oncoPredict)
## Figures: data source for Fig 3D, Fig 3E, and Fig S3C
## Dataset: GSE138709 (5 tumor + 3 adjacent iCCA samples)
## ============================================================
## Prerequisites:
##   - scRNA1: Seurat object with article-standard curate_v1 annotation
##   - GDSC2_Expr_short.rds and GDSC2_Res.rds
## Output:
##   - scRNA1 with Cisplatin_pred_group / Gemcitabine_pred_group metadata
##   - scRNA1 with lnIC50_* and IC50uM_* metadata columns
## Notes:
##   - Sensitive/resistant groups are assigned by the median predicted lnIC50,
##     preserving the supplied analysis logic.
##   - PDF Fig 3D reports the resulting cisplatin lnIC50 cutoff as 3.365177.
##   - PDF Fig 3E reports the resulting gemcitabine lnIC50 cutoff as -1.230213.
##     These values are audit targets, not hard-coded thresholds.
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(oncoPredict)
  library(dplyr)
  library(ggplot2)
})

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
set.seed(1234)

stopifnot(exists("scRNA1"), inherits(scRNA1, "Seurat"))

# ============================================================
# From Script 6.R: Load GDSC2 training data
# ============================================================
# Adjust paths to where you have the GDSC2 data files
gdsc2_expr <- readRDS("GDSC2_Expr_short.rds")   # training expression (gene x sample)
gdsc2_res  <- readRDS("GDSC2_Res.rds")           # training IC50 matrix

# ============================================================
# From Script 6.R: Extract test expression from Seurat object
# ============================================================
# Use aneuploid cells (group_copykat) as the test set
if ("group_copykat" %in% colnames(scRNA1@meta.data)) {
  cells_use <- rownames(scRNA1@meta.data)[
    !is.na(scRNA1$group_copykat) & scRNA1$group_copykat %in% c("aneuploid","diploid")
  ]
} else {
  cells_use <- colnames(scRNA1)
}

DefaultAssay(scRNA1) <- "RNA"
testExprData <- GetAssayData(scRNA1, slot = "data")[, cells_use, drop = FALSE]

# Match genes between training and test
common_genes <- intersect(rownames(gdsc2_expr), rownames(testExprData))
trainingExprData2 <- gdsc2_expr[common_genes, , drop = FALSE]
testExprData2     <- testExprData[common_genes, , drop = FALSE]

# Match drug phenotype columns
trainingPtype <- gdsc2_res

# ============================================================
# From Script 6.R: Run calcPhenotype (parameters unchanged)
# ============================================================
if (!dir.exists("calcPhenotype_Output")) dir.create("calcPhenotype_Output", recursive = TRUE)

calcPhenotype(
  trainingExprData        = as.matrix(trainingExprData2),
  trainingPtype           = as.matrix(trainingPtype),
  testExprData            = as.matrix(testExprData2),
  batchCorrect            = "standardize",
  powerTransformPhenotype = TRUE,
  removeLowVaryingGenes   = 0.2,
  minNumSamples           = 10,
  printOutput             = TRUE,
  outDir                  = "calcPhenotype_Output"
)

message("calcPhenotype complete. Results in: calcPhenotype_Output/")

# ============================================================
# From Script 7.R / Script 8.R: Read predictions and write to Seurat metadata
# ============================================================
pred <- read.csv(
  file.path("calcPhenotype_Output", "DrugPredictions.csv"),
  row.names = 1, check.names = FALSE, stringsAsFactors = FALSE
)

# Align to cells_use (keep only matching barcodes)
keep_cells <- intersect(cells_use, rownames(pred))
pred <- pred[keep_cells, , drop = FALSE]

# Build list-column: each cell gets a named numeric vector (drug -> IC50 uM)
cells_all <- rownames(scRNA1@meta.data)
vec_list  <- vector("list", length(cells_all))
names(vec_list) <- cells_all

vec_list[keep_cells] <- lapply(keep_cells, function(cb) {
  x <- as.numeric(pred[cb, ])
  names(x) <- colnames(pred)
  x
})

scRNA1@meta.data$DrugPredictions <- vec_list       # list-column in meta
scRNA1@misc$DrugPredictions      <- pred           # full matrix in misc

# ============================================================
# From Script 7.R / Script 8.R: Helper to extract drug IC50 as flat column
# ============================================================
pull_drug_to_meta <- function(obj, drug, log_scale = TRUE, prefix = NULL) {
  stopifnot("DrugPredictions" %in% colnames(obj@meta.data))
  if (is.null(prefix)) prefix <- if (log_scale) "lnIC50_" else "IC50uM_"
  vals <- vapply(obj@meta.data$DrugPredictions, function(x) {
    if (length(x) && drug %in% names(x)) x[[drug]] else NA_real_
  }, numeric(1))
  if (log_scale) vals <- log(pmax(vals, 1e-9))
  colname <- paste0(prefix, make.names(drug))
  obj[[colname]] <- vals
  message("Added meta column: ", colname)
  obj
}

# Extract Cisplatin lnIC50
# Identify the correct Cisplatin column name from predictions
cisplatin_key <- grep("cisplatin|cddp", colnames(pred), ignore.case = TRUE, value = TRUE)[1]
if (!is.na(cisplatin_key)) {
  scRNA1 <- pull_drug_to_meta(scRNA1, drug = cisplatin_key, log_scale = TRUE,  prefix = "lnIC50_")
  scRNA1 <- pull_drug_to_meta(scRNA1, drug = cisplatin_key, log_scale = FALSE, prefix = "IC50uM_")
  message("Cisplatin column used: ", cisplatin_key)
} else {
  warning("Cisplatin column not found in drug predictions.")
}

# Extract Gemcitabine lnIC50
gemcitabine_key <- grep("gemcitabine", colnames(pred), ignore.case = TRUE, value = TRUE)[1]
if (!is.na(gemcitabine_key)) {
  scRNA1 <- pull_drug_to_meta(scRNA1, drug = gemcitabine_key, log_scale = TRUE,  prefix = "lnIC50_")
  scRNA1 <- pull_drug_to_meta(scRNA1, drug = gemcitabine_key, log_scale = FALSE, prefix = "IC50uM_")
  message("Gemcitabine column used: ", gemcitabine_key)
}

# ============================================================
# From Script 7.R: Assign Cisplatin sensitivity groups (aneuploid cells)
# ============================================================
# Find the lnIC50_Cisplatin column
fetch_ic50_for_drug <- function(obj, drug_pattern) {
  mdn <- colnames(obj@meta.data)
  # 1) direct ln column
  hit_ln <- grep(paste0("^lnIC50.*", drug_pattern), mdn, ignore.case = TRUE, value = TRUE)
  if (length(hit_ln)) return(list(col = hit_ln[1], vals = obj@meta.data[[hit_ln[1]]]))
  # 2) IC50uM -> ln
  hit_raw <- grep(paste0("^IC50uM.*", drug_pattern), mdn, ignore.case = TRUE, value = TRUE)
  if (length(hit_raw)) {
    v <- log(pmax(as.numeric(obj@meta.data[[hit_raw[1]]]), 1e-9))
    return(list(col = hit_raw[1], vals = v))
  }
  # 3) DrugPredictions list-column
  if ("DrugPredictions" %in% mdn) {
    all_keys <- unique(unlist(lapply(obj@meta.data$DrugPredictions, names)))
    k <- grep(drug_pattern, all_keys, ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(k)) {
      vals <- vapply(obj@meta.data$DrugPredictions, function(x) {
        if (length(x) && k %in% names(x)) as.numeric(x[[k]]) else NA_real_
      }, numeric(1))
      return(list(col = k, vals = log(pmax(vals, 1e-9))))
    }
  }
  # 4) misc table
  if (!is.null(obj@misc$DrugPredictions)) {
    dp <- obj@misc$DrugPredictions
    k  <- grep(drug_pattern, colnames(dp), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(k)) {
      v <- rep(NA_real_, ncol(obj)); names(v) <- colnames(obj)
      common <- intersect(names(v), rownames(dp))
      v[common] <- as.numeric(dp[common, k])
      return(list(col = k, vals = log(pmax(v, 1e-9))))
    }
  }
  # 5) CSV fallback
  fcsv <- file.path("calcPhenotype_Output", "DrugPredictions.csv")
  if (file.exists(fcsv)) {
    dp <- read.csv(fcsv, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
    k  <- grep(drug_pattern, colnames(dp), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(k)) {
      v <- rep(NA_real_, ncol(obj)); names(v) <- colnames(obj)
      common <- intersect(names(v), rownames(dp))
      v[common] <- as.numeric(dp[common, k])
      return(list(col = k, vals = log(pmax(v, 1e-9))))
    }
  }
  stop("Cannot find IC50 for drug pattern: ", drug_pattern)
}

# Assign Cisplatin_pred_group
assign_drug_group <- function(obj, drug_pattern, group_colname, subset_cells = NULL) {
  ic50_info <- fetch_ic50_for_drug(obj, drug_pattern)
  ic50_vals <- ic50_info$vals
  names(ic50_vals) <- rownames(obj@meta.data)

  if (!is.null(subset_cells)) {
    ic50_vals_sub <- ic50_vals[subset_cells]
    thr <- median(ic50_vals_sub, na.rm = TRUE)
  } else {
    thr <- median(ic50_vals, na.rm = TRUE)
  }

  grp <- ifelse(ic50_vals >= thr, "Predicted resistant", "Predicted sensitive")
  grp[!is.finite(ic50_vals)] <- NA_character_
  obj[[group_colname]] <- factor(grp, levels = c("Predicted sensitive", "Predicted resistant"))
  message("Assigned ", group_colname, " (median threshold: ", round(thr, 4), ")")
  if (grepl("cisplatin|cddp", drug_pattern, ignore.case = TRUE)) {
    message("  PDF Fig 3D reported lnIC50 cutoff: 3.365177 (audit target).")
  }
  if (grepl("gemcitabine", drug_pattern, ignore.case = TRUE)) {
    message("  PDF Fig 3E reported lnIC50 cutoff: -1.230213 (audit target).")
  }
  obj
}

# Assign Cisplatin group (all cells)
scRNA1 <- assign_drug_group(scRNA1, "cisplatin|cddp", "Cisplatin_pred_group")

# Assign Gemcitabine group (all cells)
scRNA1 <- assign_drug_group(scRNA1, "gemcitabine", "Gemcitabine_pred_group")

# ============================================================
# From Script 5.R: QC VlnPlot for IC50
# ============================================================
if ("lnIC50_Cisplatin" %in% colnames(scRNA1@meta.data) ||
    any(grepl("lnIC50.*Cisplatin", colnames(scRNA1@meta.data)))) {
  ln_col <- grep("lnIC50.*Cisplatin", colnames(scRNA1@meta.data), value = TRUE)[1]
  p_vln_ic50 <- VlnPlot(
    scRNA1,
    features = ln_col,
    group.by = "Cisplatin_pred_group",
    pt.size  = 0
  )
  if (!dir.exists("output")) dir.create("output", recursive = TRUE)
  ggsave("output/VlnPlot_lnIC50_Cisplatin_by_pred_group.pdf", p_vln_ic50, width = 6, height = 4)
}

# ============================================================
# From Script 7.R: Save object with drug predictions
# ============================================================
saveRDS(scRNA1, file = "output/scRNA1_with_DrugPredictions.rds")

message("Done: 03_drug_prediction.R complete.")
message("  Drug group columns added: Cisplatin_pred_group, Gemcitabine_pred_group")
