## ============================================================
## YBX1 Regulon AUC Density UMAP
## Figures: Fig S3D
## Dataset: GSE138709 (5 tumor + 3 adjacent iCCA samples)
## ============================================================
## Prerequisites:
##   - scRNAauc or scRNA1 with YBX1 regulon AUC metadata from
##     09_data_SCENIC_AUC_integration_for_Fig3_FigS3.R
## Output:
##   - output/SCENIC_FeatureDensity_<YBX1_regulon_AUC>_umap.pdf/.png (Fig S3D)
## Notes:
##   - This script is restricted to Fig S3D only. The cisplatin heatmap belongs
##     to Fig S3B, and the cisplatin IC50 density belongs to Fig S3C.
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Nebulosa)
  library(Matrix)
  library(ggplot2)
})

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

if (exists("scRNAauc") && inherits(scRNAauc, "Seurat")) {
  obj <- scRNAauc
} else if (exists("scRNA1") && inherits(scRNA1, "Seurat")) {
  obj <- scRNA1
} else {
  candidates <- c(
    "output/scRNA1_with_SCENIC_AUC.rds",
    "output/scRNA1_with_DrugPredictions.rds",
    "output/scRNA1_annotated.rds",
    "output/scRNA1_preprocessed.rds"
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) {
    stop("No Seurat object was found in memory and no expected RDS file exists in output/.")
  }
  obj <- readRDS(hit)
}

dir.create("output", showWarnings = FALSE, recursive = TRUE)

tf_root <- function(x) {
  y <- sub("\\s*\\(\\d+g\\)$", "", x)
  y <- sub("_[0-9]+g$", "", y, ignore.case = TRUE)
  y <- sub("_extended$", "", y, ignore.case = TRUE)
  y
}

load_auc_matrix <- function() {
  if (exists("auc_mtx")) return(auc_mtx)
  cand <- list.files("int", pattern = "^3\\.4_regulonAUC.*\\.Rds$",
                     full.names = TRUE)
  if (!length(cand)) return(NULL)
  if (!requireNamespace("SCENIC", quietly = TRUE)) {
    warning("SCENIC is not installed, so regulonAUC RDS cannot be converted with getAUC().")
    return(NULL)
  }
  SCENIC::getAUC(readRDS(cand[1]))
}

find_ybx1_regulon <- function(obj) {
  mdn <- colnames(obj@meta.data)

  preferred <- intersect(c("YBX1_extended_907g", "YBX1_extended_907G"), mdn)
  if (length(preferred)) {
    return(list(obj = obj, col = preferred[1]))
  }

  md_hits <- grep("^YBX1.*(regulon|AUC|extended|[0-9]+g)", mdn,
                  ignore.case = TRUE, value = TRUE)
  md_hits <- md_hits[tf_root(md_hits) == "YBX1"]
  if (length(md_hits)) {
    non_ext <- md_hits[!grepl("_extended$", md_hits, ignore.case = TRUE)]
    return(list(obj = obj, col = if (length(non_ext)) non_ext[1] else md_hits[1]))
  }

  auc <- load_auc_matrix()
  if (!is.null(auc)) {
    common <- intersect(colnames(obj), colnames(auc))
    if (!length(common)) stop("YBX1 regulon AUC matrix and Seurat object share no cells.")

    auc_hits <- rownames(auc)[tf_root(rownames(auc)) == "YBX1"]
    if (!length(auc_hits)) stop("YBX1 regulon was not found in auc_mtx rownames.")

    preferred_auc <- auc_hits[auc_hits == "YBX1_extended_907g"]
    pick <- if (length(preferred_auc)) preferred_auc[1] else {
      non_ext <- auc_hits[!grepl("_extended$", auc_hits, ignore.case = TRUE)]
      if (length(non_ext)) non_ext[1] else auc_hits[1]
    }

    vals <- rep(NA_real_, ncol(obj))
    names(vals) <- colnames(obj)
    vals[common] <- as.numeric(auc[pick, common])
    col <- make.names(pick)
    obj[[col]] <- vals
    return(list(obj = obj, col = col))
  }

  stop("No YBX1 regulon AUC metadata or AUC matrix was found.")
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

ybx1 <- find_ybx1_regulon(obj)
obj <- ybx1$obj
ybx1_col <- ybx1$col
message("Fig S3D YBX1 regulon AUC feature: ", ybx1_col)

red_use <- if ("umap" %in% names(obj@reductions)) {
  "umap"
} else if ("tsne" %in% names(obj@reductions)) {
  "tsne"
} else {
  stop("The Seurat object has neither UMAP nor tSNE reduction.")
}

assay_result <- metadata_feature_to_assay(obj, ybx1_col, "YBX1_Regulon_AUC")
obj <- assay_result$obj
plot_feature <- assay_result$feature

p_s3d <- Nebulosa::plot_density(
  obj,
  features = plot_feature,
  reduction = red_use
) + ggtitle("YB-1 regulon AUC")

ggsave(file.path("output", paste0("SCENIC_FeatureDensity_", plot_feature, "_", red_use, ".pdf")),
       p_s3d, width = 7, height = 6)
ggsave(file.path("output", paste0("SCENIC_FeatureDensity_", plot_feature, "_", red_use, ".png")),
       p_s3d, width = 7, height = 6, dpi = 300, bg = "white")

saveRDS(obj, file = "output/scRNA_with_YBX1_regulon_density_input.rds")
message("Done: Fig S3D YBX1 regulon AUC density UMAP saved.")
