## ============================================================
## YBX1 Expression, Regulon AUC, and GSVA UMAP Visualizations
## Figures: Fig 1C (YB-1 density UMAP), Fig 1F (YB-1 expression UMAP),
##          Fig 3B top (YB-1 regulon AUC UMAP),
##          Fig 3C top (GSVA YB-1 target score UMAP)
## Dataset: GSE138709 (5 tumor + 3 adjacent iCCA samples)
## ============================================================
## Prerequisites:
##   - scRNA1 with group_copykat metadata
##     produced by 09_data_SCENIC_AUC_integration_for_Fig3_FigS3.R
## Output:
##   - FeatureDensity_YBX1_umap.pdf/.png            (Fig 1C)
##   - FeaturePlot_YBX1_umap.pdf/.png               (Fig 1F)
##   - FeaturePlot_YBX1_extended_907g_umap.pdf/.png (Fig 3B top)
##   - FeaturePlot_GSVA_YBX1_targets_ssgsea_umap.pdf/.png (Fig 3C top)
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

## Pre-flight checks
stopifnot(exists("scRNA1"), inherits(scRNA1, "Seurat"))
stopifnot("GSVA_YBX1_targets_ssgsea" %in% colnames(scRNA1@meta.data))
stopifnot("YBX1_extended_907g" %in% colnames(scRNA1@meta.data))

## 若默认 assay 没有 YBX1，则切换到含有 YBX1 的 assay（不改作图参数）
has_ybx1 <- function(obj, assay) {
  "YBX1" %in% tryCatch(rownames(obj@assays[[assay]]), error = function(e) character(0))
}
if (!has_ybx1(scRNA1, DefaultAssay(scRNA1))) {
  a_found <- NA_character_
  for (a in names(scRNA1@assays)) if (has_ybx1(scRNA1, a)) { a_found <- a; break }
  if (is.na(a_found)) stop("未在任何 assay 中找到基因 YBX1。")
  DefaultAssay(scRNA1) <- a_found
}

## Fig 1C: density-style YB-1 UMAP, matching the PDF panel.
if (!requireNamespace("Nebulosa", quietly = TRUE)) {
  stop("Fig 1C requires the Nebulosa package for density UMAP plotting.")
}
p0 <- Nebulosa::plot_density(
  scRNA1,
  features = "YBX1",
  reduction = "umap"
) + ggtitle("YB-1")

## Feature plots for Fig 1F, Fig 3B top, and Fig 3C top.
p1 <- FeaturePlot(scRNA1, features = "YBX1_extended_907g", label = F, reduction = "umap")
p2 <- FeaturePlot(scRNA1, features = "GSVA_YBX1_targets_ssgsea", label = F, reduction = "umap")
p3 <- FeaturePlot(scRNA1, features = "YBX1", label = F, reduction = "umap")

## 单独保存
if (!dir.exists("output")) dir.create("output", recursive = TRUE)

ggsave("output/FeatureDensity_YBX1_umap.pdf", p0, width = 7, height = 6)
ggsave("output/FeatureDensity_YBX1_umap.png", p0, width = 7, height = 6, dpi = 300, bg = "white")

ggsave("output/FeaturePlot_YBX1_extended_907g_umap.pdf", p1, width = 7, height = 6)
ggsave("output/FeaturePlot_YBX1_extended_907g_umap.png", p1, width = 7, height = 6, dpi = 300, bg = "white")

ggsave("output/FeaturePlot_GSVA_YBX1_targets_ssgsea_umap.pdf", p2, width = 7, height = 6)
ggsave("output/FeaturePlot_GSVA_YBX1_targets_ssgsea_umap.png", p2, width = 7, height = 6, dpi = 300, bg = "white")

ggsave("output/FeaturePlot_YBX1_umap.pdf", p3, width = 7, height = 6)
ggsave("output/FeaturePlot_YBX1_umap.png", p3, width = 8, height = 6, dpi = 300, bg = "white")

# 如需在会话里显示：
p0; p1; p2; p3
