## ============================================================
## SCENIC Regulon Volcano Plot: Cisplatin Sensitive vs Resistant
## Figures: Fig S3A
## Dataset: GSE138709 (5 tumor + 3 adjacent iCCA samples)
## ============================================================
## Prerequisites:
##   - auc_mtx or output/auc_mtx_regulon_by_cell.rds
##   - scRNAauc/scRNA1 with Cisplatin_pred_group metadata
## Output:
##   - output/SCENIC_RegulonAUC_Diff_Cisplatin_pred_group.csv
##   - output/SCENIC_Regulon_Volcano_Cisplatin_pred_group.pdf/.png
## Notes:
##   - Fig S3A and Fig S3B use the same Cisplatin_pred_group.
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
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
  stop("Both Predicted sensitive and Predicted resistant cells are required for Fig S3A.")
}

mean_s <- rowMeans(auc_mtx[, idx_s, drop = FALSE], na.rm = TRUE)
mean_r <- rowMeans(auc_mtx[, idx_r, drop = FALSE], na.rm = TRUE)
mean_diff <- mean_r - mean_s

pvals <- apply(auc_mtx, 1, function(v) {
  tryCatch(
    stats::wilcox.test(v[idx_r], v[idx_s], exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
})
fdr <- p.adjust(pvals, method = "BH")

df <- data.frame(
  Regulon = rownames(auc_mtx),
  mean_sensitive = mean_s,
  mean_resistant = mean_r,
  meanDiff = mean_diff,
  FDR = fdr,
  check.names = FALSE
)
df$TF_root <- tf_root(df$Regulon)
df$is_extended <- grepl("_extended$", df$Regulon, ignore.case = TRUE)
df <- df[order(is.na(df$FDR), df$FDR, -abs(df$meanDiff), df$is_extended), ]
df_dedup <- df[!duplicated(df$TF_root), , drop = FALSE]

write.csv(
  df_dedup,
  file.path("output", paste0("SCENIC_RegulonAUC_Diff_Cisplatin_pred_group", suffix, ".csv")),
  row.names = FALSE
)

df_dedup$direction <- "Not significant"
df_dedup$direction[df_dedup$FDR < 0.05 & df_dedup$meanDiff > 0.001] <- "Resistant"
df_dedup$direction[df_dedup$FDR < 0.05 & df_dedup$meanDiff < -0.001] <- "Sensitive"
df_dedup$negLogFDR <- -log10(pmax(df_dedup$FDR, .Machine$double.xmin))

top_labels <- head(df_dedup$Regulon[order(df_dedup$FDR, -abs(df_dedup$meanDiff))], 10)
ybx1_label <- df_dedup$Regulon[df_dedup$TF_root == "YBX1"][1]
label_names <- unique(na.omit(c(top_labels, ybx1_label)))
label_df <- df_dedup[df_dedup$Regulon %in% label_names, , drop = FALSE]

p <- ggplot(df_dedup, aes(x = meanDiff, y = negLogFDR)) +
  geom_point(aes(color = direction), alpha = 0.7, size = 2) +
  geom_vline(xintercept = c(-0.001, 0.001), color = "grey40",
             linetype = "longdash", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), color = "grey40",
             linetype = "longdash", linewidth = 0.4) +
  geom_point(data = label_df, shape = 1, color = "black", size = 4, stroke = 0.8) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = Regulon),
    size = 4,
    max.overlaps = 1000,
    box.padding = 0.35,
    point.padding = 0.3
  ) +
  scale_color_manual(values = c(
    "Sensitive" = "#329E3F",
    "Resistant" = "#ED4F4F",
    "Not significant" = "grey70"
  )) +
  labs(
    x = "AUC mean difference: resistant - sensitive",
    y = "-log10(FDR)",
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank())

ggsave(file.path("output", paste0("SCENIC_Regulon_Volcano_Cisplatin_pred_group", suffix, ".pdf")),
       p, width = 7, height = 6)
ggsave(file.path("output", paste0("SCENIC_Regulon_Volcano_Cisplatin_pred_group", suffix, ".png")),
       p, width = 7, height = 6, dpi = 300, bg = "white")

message("Done: Fig S3A cisplatin regulon volcano saved.")
