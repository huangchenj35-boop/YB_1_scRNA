## ============================================================
## Cisplatin IC50 UMAP Density Plot
## Figures: Fig S3C
## Dataset: GSE138709 (5 tumor + 3 adjacent iCCA samples)
## ============================================================
## Prerequisites:
##   - scRNA1 with cisplatin prediction metadata from
##     07_data_drug_sensitivity_prediction_for_Fig3D_Fig3E_FigS3C.R
## Output:
##   - output/FeatureDensity_<cisplatin_lnIC50>_umap.pdf/.png (Fig S3C)
## Notes:
##   - Sensitive/resistant grouping, when recreated here for audit, uses the
##     original median lnIC50 rule. The PDF Fig 3D cutoff value 3.365177 is an
##     audit target only and is not hard-coded as the grouping threshold.
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Nebulosa)
  library(Matrix)
  library(ggplot2)
})

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

if (!exists("scRNA1") || !inherits(scRNA1, "Seurat")) {
  candidates <- c(
    "output/scRNA1_with_DrugPredictions.rds",
    "output/scRNA1_annotated.rds",
    "output/scRNA1_preprocessed.rds"
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) {
    stop("scRNA1 was not found in memory and no expected RDS file exists in output/.")
  }
  scRNA1 <- readRDS(hit)
}

dir.create("output", showWarnings = FALSE, recursive = TRUE)

find_cisplatin_ln_ic50 <- function(obj) {
  mdn <- colnames(obj@meta.data)

  ln_hit <- grep("^(lnIC50|logIC50).*?(cisplatin|cddp)", mdn,
                 ignore.case = TRUE, value = TRUE)
  if (length(ln_hit)) {
    return(list(obj = obj, col = ln_hit[1], values = obj@meta.data[[ln_hit[1]]]))
  }

  raw_hit <- grep("^(IC50uM|IC50).*?(cisplatin|cddp)", mdn,
                  ignore.case = TRUE, value = TRUE)
  raw_hit <- raw_hit[!grepl("^(lnIC50|logIC50)", raw_hit, ignore.case = TRUE)]
  if (length(raw_hit)) {
    new_col <- paste0("lnIC50_", sub("^(IC50uM_|IC50_)", "", raw_hit[1]))
    obj[[new_col]] <- log(pmax(as.numeric(obj@meta.data[[raw_hit[1]]]), 1e-9))
    return(list(obj = obj, col = new_col, values = obj@meta.data[[new_col]]))
  }

  if (!is.null(obj@misc$DrugPredictions)) {
    pred <- obj@misc$DrugPredictions
  } else if (file.exists(file.path("calcPhenotype_Output", "DrugPredictions.csv"))) {
    pred <- read.csv(file.path("calcPhenotype_Output", "DrugPredictions.csv"),
                     row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    pred <- NULL
  }

  if (!is.null(pred)) {
    drug_col <- grep("cisplatin|cddp", colnames(pred), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(drug_col)) {
      vals_raw <- rep(NA_real_, ncol(obj))
      names(vals_raw) <- colnames(obj)
      common <- intersect(names(vals_raw), rownames(pred))
      vals_raw[common] <- as.numeric(pred[common, drug_col])
      new_col <- paste0("lnIC50_", make.names(drug_col))
      obj[[new_col]] <- log(pmax(vals_raw, 1e-9))
      return(list(obj = obj, col = new_col, values = obj@meta.data[[new_col]]))
    }
  }

  stop("No cisplatin lnIC50 source was found in metadata, obj@misc, or calcPhenotype_Output/DrugPredictions.csv.")
}

safe_feature_name <- function(x) {
  make.names(x, unique = TRUE)
}

metadata_feature_to_assay <- function(obj, meta_col, assay_name) {
  stopifnot(meta_col %in% colnames(obj@meta.data))
  feature <- safe_feature_name(meta_col)
  vals <- as.numeric(obj@meta.data[[meta_col]])
  vals[!is.finite(vals)] <- 0
  mat <- Matrix::Matrix(vals, nrow = 1, sparse = TRUE)
  rownames(mat) <- feature
  colnames(mat) <- colnames(obj)
  obj[[assay_name]] <- CreateAssayObject(data = mat)
  DefaultAssay(obj) <- assay_name
  list(obj = obj, feature = feature)
}

ic50 <- find_cisplatin_ln_ic50(scRNA1)
scRNA1 <- ic50$obj
ic50_col <- ic50$col
ic50_vals <- as.numeric(scRNA1@meta.data[[ic50_col]])

if (!"Cisplatin_pred_group" %in% colnames(scRNA1@meta.data)) {
  thr <- stats::median(ic50_vals, na.rm = TRUE)
  if (!is.finite(thr)) stop("Median cisplatin lnIC50 is not finite.")
  grp <- ifelse(ic50_vals >= thr, "Predicted resistant", "Predicted sensitive")
  grp[!is.finite(ic50_vals)] <- NA_character_
  scRNA1[["Cisplatin_pred_group"]] <- factor(
    grp,
    levels = c("Predicted sensitive", "Predicted resistant")
  )
  message("Cisplatin_pred_group recreated using median lnIC50 threshold: ", signif(thr, 6))
}

message("PDF Fig 3D reported cisplatin lnIC50 cutoff for audit: 3.365177")
message("Fig S3C density feature: ", ic50_col)
message("Cisplatin_pred_group counts:")
print(table(scRNA1$Cisplatin_pred_group, useNA = "ifany"))

red_use <- if ("umap" %in% names(scRNA1@reductions)) {
  "umap"
} else if ("tsne" %in% names(scRNA1@reductions)) {
  "tsne"
} else {
  stop("scRNA1 has neither UMAP nor tSNE reduction.")
}

assay_result <- metadata_feature_to_assay(scRNA1, ic50_col, "Cisplatin_lnIC50")
scRNA1 <- assay_result$obj
plot_feature <- assay_result$feature

p_s3c <- Nebulosa::plot_density(
  scRNA1,
  features = plot_feature,
  reduction = red_use
) + ggtitle("Cisplatin lnIC50")

ggsave(file.path("output", paste0("FeatureDensity_", plot_feature, "_", red_use, ".pdf")),
       p_s3c, width = 7, height = 6)
ggsave(file.path("output", paste0("FeatureDensity_", plot_feature, "_", red_use, ".png")),
       p_s3c, width = 7, height = 6, dpi = 300, bg = "white")

saveRDS(scRNA1, file = "output/scRNA1_with_cisplatin_lnIC50_density_input.rds")
message("Done: Fig S3C cisplatin lnIC50 density UMAP saved.")
