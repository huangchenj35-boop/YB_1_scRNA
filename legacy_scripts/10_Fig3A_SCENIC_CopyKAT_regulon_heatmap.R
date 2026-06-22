## ============================================================
## SCENIC Regulon Heatmap: Aneuploid vs Diploid
## Figures: Fig 3A
## Dataset: GSE138709 (5 tumor + 3 adjacent iCCA samples)
## ============================================================
## Prerequisites:
##   - auc_mtx or output/auc_mtx_regulon_by_cell.rds
##   - scRNAauc/scRNA1 with group_copykat metadata
## Output:
##   - output/SCENIC_RegulonAUC_Diff_CopyKat2groups.csv
##   - output/SCENIC_RegulonAUC_Top50_GroupMean_heatmap_CopyKat2groups.pdf/.png
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

if (exists("scRNAauc") && inherits(scRNAauc, "Seurat")) {
  scRNA <- scRNAauc
} else if (exists("scRNA1") && inherits(scRNA1, "Seurat")) {
  scRNA <- scRNA1
} else if (file.exists("output/scRNAauc.rds")) {
  scRNA <- readRDS("output/scRNAauc.rds")
} else {
  stop("No scRNAauc/scRNA1 object was found.")
}

if (!"group_copykat" %in% colnames(scRNA@meta.data)) {
  stop("group_copykat metadata is missing. Run 04_Fig1E_FigS2_inferCNV_CopyKAT.R first.")
}

tf_root <- function(x) {
  y <- sub("\\s*\\(\\d+g\\)$", "", x)
  y <- sub("_[0-9]+g$", "", y, ignore.case = TRUE)
  y <- sub("_extended$", "", y, ignore.case = TRUE)
  y
}

dedup_regulons <- function(df) {
  df$TF_root <- tf_root(df$regulon)
  df <- df[order(is.na(df$FDR), df$FDR, df$p.value, -abs(df$diff)), ]
  df[!duplicated(df$TF_root), , drop = FALSE]
}

common <- intersect(colnames(auc_mtx), colnames(scRNA))
if (!length(common)) stop("auc_mtx and Seurat object share no cells.")
auc_mtx <- auc_mtx[, common, drop = FALSE]

group <- factor(scRNA@meta.data[common, "group_copykat"],
                levels = c("aneuploid", "diploid"))
keep <- !is.na(group)
auc_mtx <- auc_mtx[, keep, drop = FALSE]
group <- droplevels(group[keep])

idx_a <- which(group == "aneuploid")
idx_d <- which(group == "diploid")
if (!length(idx_a) || !length(idx_d)) {
  stop("Both aneuploid and diploid cells are required for Fig 3A.")
}

mean_a <- rowMeans(auc_mtx[, idx_a, drop = FALSE], na.rm = TRUE)
mean_d <- rowMeans(auc_mtx[, idx_d, drop = FALSE], na.rm = TRUE)
diff_v <- mean_a - mean_d

pvals <- apply(auc_mtx, 1, function(v) {
  tryCatch(
    stats::wilcox.test(v[idx_a], v[idx_d], exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
})
fdr <- p.adjust(pvals, method = "BH")

diff_res <- data.frame(
  regulon = rownames(auc_mtx),
  mean_aneuploid = mean_a,
  mean_diploid = mean_d,
  diff = diff_v,
  p.value = pvals,
  FDR = fdr,
  check.names = FALSE
)
diff_res <- diff_res[order(is.na(diff_res$FDR), diff_res$FDR, -abs(diff_res$diff)), ]

write.csv(
  diff_res,
  file.path("output", paste0("SCENIC_RegulonAUC_Diff_CopyKat2groups", suffix, ".csv")),
  row.names = FALSE
)

dedup <- dedup_regulons(diff_res)
top <- subset(dedup, diff > 0 & is.finite(FDR) & is.finite(p.value))
top <- top[order(top$FDR, top$p.value, -top$diff), , drop = FALSE]
top50 <- head(top$regulon, 50)
if (length(top50) < 2) stop("Fewer than two up-in-aneuploid regulons were available for Fig 3A.")

ht_mat <- cbind(
  aneuploid = rowMeans(auc_mtx[top50, idx_a, drop = FALSE], na.rm = TRUE),
  diploid = rowMeans(auc_mtx[top50, idx_d, drop = FALSE], na.rm = TRUE)
)
rownames(ht_mat) <- top50

anno_col <- data.frame(
  group_copykat = factor(c("aneuploid", "diploid"), levels = c("aneuploid", "diploid"))
)
rownames(anno_col) <- colnames(ht_mat)

bk <- c(seq(-3, -0.1, by = 0.01), seq(0, 3, by = 0.01))
heat_colors <- c(
  colorRampPalette(c("#2166ac", "#f7fbff"))(length(bk) / 2),
  colorRampPalette(c("#f7fbff", "#b2182b"))(length(bk) / 2)
)
annotation_colors <- list(
  group_copykat = c("aneuploid" = "#F29C64", "diploid" = "#5AAFD6")
)

pdf(file.path("output", paste0("SCENIC_RegulonAUC_Top50_GroupMean_heatmap_CopyKat2groups", suffix, ".pdf")),
    width = 8, height = 12)
pheatmap::pheatmap(
  ht_mat,
  annotation_col = anno_col,
  annotation_colors = annotation_colors,
  show_colnames = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  color = heat_colors,
  breaks = bk,
  scale = "row",
  legend_breaks = seq(-2, 2, 2),
  annotation_names_row = FALSE,
  annotation_names_col = TRUE,
  fontsize_row = 8
)
dev.off()

png(file.path("output", paste0("SCENIC_RegulonAUC_Top50_GroupMean_heatmap_CopyKat2groups", suffix, ".png")),
    width = 2400, height = 3600, res = 300)
pheatmap::pheatmap(
  ht_mat,
  annotation_col = anno_col,
  annotation_colors = annotation_colors,
  show_colnames = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  color = heat_colors,
  breaks = bk,
  scale = "row",
  legend_breaks = seq(-2, 2, 2),
  annotation_names_row = FALSE,
  annotation_names_col = TRUE,
  fontsize_row = 8
)
dev.off()

message("Done: Fig 3A CopyKAT regulon heatmap saved.")
