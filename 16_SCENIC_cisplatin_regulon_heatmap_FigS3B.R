## ============================================================
## SCENIC Regulon AUC Heatmap: Cisplatin Sensitive vs Resistant
## Figures: Fig S3B
## Dataset: GSE138709 (5 tumor + 3 adjacent iCCA samples)
## ============================================================
## Panel note:
##   Fig S3A and Fig S3B use the same grouping:
##   Cisplatin_pred_group = Predicted sensitive vs Predicted resistant.
##   Fig S3A is the regulon volcano plot; Fig S3B is the matching
##   regulon AUC heatmap.
## Output:
##   - output/SCENIC_RegulonAUC_Top20_GroupMean_heatmap_CisplatinGroup_Cisplatin.pdf/.png
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(pheatmap)
  library(RColorBrewer)
})

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
if (!exists("suffix")) suffix <- ""
dir.create("output", showWarnings = FALSE, recursive = TRUE)

if (!exists("auc_mtx")) {
  if (file.exists("output/auc_mtx_regulon_by_cell.rds")) {
    auc_mtx <- readRDS("output/auc_mtx_regulon_by_cell.rds")
  } else {
    stop("auc_mtx was not found. Run 09_data_SCENIC_AUC_integration_for_Fig3_FigS3.R first.")
  }
}

if (exists("scRNAauc_mc") && inherits(scRNAauc_mc, "Seurat")) {
  scRNA <- scRNAauc_mc
} else if (exists("scRNAauc") && inherits(scRNAauc, "Seurat")) {
  scRNA <- scRNAauc
} else if (exists("scRNA1") && inherits(scRNA1, "Seurat")) {
  scRNA <- scRNA1
} else if (file.exists("output/scRNAauc.rds")) {
  scRNA <- readRDS("output/scRNAauc.rds")
} else {
  stop("No scRNAauc_mc/scRNAauc/scRNA1 object was found.")
}

if (!"Cisplatin_pred_group" %in% colnames(scRNA@meta.data)) {
  stop("Cisplatin_pred_group metadata is missing. Run 07_data_drug_sensitivity_prediction_for_Fig3D_Fig3E_FigS3C.R first.")
}

tf_root <- function(x) {
  y <- sub("\\s*\\(\\d+g\\)$", "", x)
  y <- sub("_[0-9]+g$", "", y, ignore.case = TRUE)
  y <- sub("_extended$", "", y, ignore.case = TRUE)
  y
}

common <- intersect(colnames(auc_mtx), colnames(scRNA))
if (!length(common)) stop("auc_mtx and Seurat object share no cells.")
auc_mtx <- auc_mtx[, common, drop = FALSE]

group <- factor(scRNA@meta.data[common, "Cisplatin_pred_group"],
                levels = c("Predicted sensitive", "Predicted resistant"))
keep <- !is.na(group)
auc_mtx <- auc_mtx[, keep, drop = FALSE]
group <- droplevels(group[keep])

idx_s <- which(group == "Predicted sensitive")
idx_r <- which(group == "Predicted resistant")
if (!length(idx_s) || !length(idx_r)) {
  stop("Both Predicted sensitive and Predicted resistant cells are required for Fig S3B.")
}

mean_s <- rowMeans(auc_mtx[, idx_s, drop = FALSE], na.rm = TRUE)
mean_r <- rowMeans(auc_mtx[, idx_r, drop = FALSE], na.rm = TRUE)
diff_v <- mean_s - mean_r

diff_res <- data.frame(
  regulon = rownames(auc_mtx),
  mean_sensitive = mean_s,
  mean_resistant = mean_r,
  diff_sensitive_minus_resistant = diff_v,
  TF_root = tf_root(rownames(auc_mtx)),
  is_extended = grepl("_extended$", rownames(auc_mtx), ignore.case = TRUE),
  check.names = FALSE
)
diff_res <- diff_res[order(-abs(diff_res$diff_sensitive_minus_resistant), diff_res$is_extended), ]
diff_res <- diff_res[!duplicated(diff_res$TF_root), , drop = FALSE]

ybx1_pick <- diff_res$regulon[diff_res$TF_root == "YBX1"][1]
top20 <- head(diff_res$regulon, 20)
if (!is.na(ybx1_pick) && !(ybx1_pick %in% top20)) {
  top20 <- unique(c(ybx1_pick, head(top20, 19)))
}
if (length(top20) < 2) stop("Fewer than two regulons were available for Fig S3B.")

write.csv(
  diff_res,
  file.path("output", paste0("SCENIC_RegulonAUC_Diff_Cisplatin_pred_group_for_heatmap", suffix, ".csv")),
  row.names = FALSE
)

ht_mat <- cbind(
  `Predicted sensitive` = rowMeans(auc_mtx[top20, idx_s, drop = FALSE], na.rm = TRUE),
  `Predicted resistant` = rowMeans(auc_mtx[top20, idx_r, drop = FALSE], na.rm = TRUE)
)
rownames(ht_mat) <- top20

anno_col <- data.frame(
  Cisplatin_pred_group = factor(
    c("Predicted sensitive", "Predicted resistant"),
    levels = c("Predicted sensitive", "Predicted resistant")
  )
)
rownames(anno_col) <- colnames(ht_mat)

bk <- c(seq(-3, -0.1, by = 0.01), seq(0, 3, by = 0.01))
heat_colors <- c(
  colorRampPalette(c("#2166ac", "#f7fbff"))(length(bk) / 2),
  colorRampPalette(c("#f7fbff", "#b2182b"))(length(bk) / 2)
)
annotation_colors <- list(
  Cisplatin_pred_group = c(
    "Predicted sensitive" = "#377EB8",
    "Predicted resistant" = "#E41A1C"
  )
)

pdf(file.path("output", paste0("SCENIC_RegulonAUC_Top20_GroupMean_heatmap_CisplatinGroup_Cisplatin", suffix, ".pdf")),
    width = 7, height = 8)
pheatmap::pheatmap(
  ht_mat,
  annotation_col = anno_col,
  annotation_colors = annotation_colors,
  show_colnames = TRUE,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "row",
  color = heat_colors,
  breaks = bk,
  border_color = NA,
  legend_breaks = seq(-2, 2, 2),
  fontsize_row = 8
)
dev.off()

png(file.path("output", paste0("SCENIC_RegulonAUC_Top20_GroupMean_heatmap_CisplatinGroup_Cisplatin", suffix, ".png")),
    width = 2100, height = 2400, res = 300)
pheatmap::pheatmap(
  ht_mat,
  annotation_col = anno_col,
  annotation_colors = annotation_colors,
  show_colnames = TRUE,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "row",
  color = heat_colors,
  breaks = bk,
  border_color = NA,
  legend_breaks = seq(-2, 2, 2),
  fontsize_row = 8
)
dev.off()

message("Done: Fig S3B cisplatin regulon heatmap saved.")
