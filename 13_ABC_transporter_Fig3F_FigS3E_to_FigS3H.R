## ============================================================
## ABC Transporter Expression and Heatmap Panels
## Figures: Fig 3F (ABCB1/ABCC1/ABCC2/MVP density UMAP),
##          Fig S3E-S3H (YB-1 + ABC transporter heatmaps)
## Dataset: GSE138709 (5 tumor + 3 adjacent iCCA samples)
## ============================================================
## Prerequisites:
##   - scRNA1 with UMAP reduction
##   - scRNA1 metadata columns:
##       sample_group, group_copykat,
##       Cisplatin_pred_group, Gemcitabine_pred_group
##   - RNA/SCT assay containing YBX1, ABCB1, ABCC1, ABCC2, and MVP
## Output:
##   - Fig3F_ABC_transporters_density_UMAP.pdf/.png
##   - FeatureDensity_ABCB1_umap.pdf/.png, etc.
##   - Heatmap_YBX1_ABC_Tumor_vs_Adjacent.pdf       (Fig S3E)
##   - Heatmap_YBX1_ABC_aneuploid_vs_diploid.pdf    (Fig S3F)
##   - Heatmap_YBX1_ABC_Cisplatin.pdf               (Fig S3G)
##   - Heatmap_YBX1_ABC_Gemcitabine.pdf             (Fig S3H)
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
})

if (!requireNamespace("Nebulosa", quietly = TRUE)) {
  stop("Fig 3F requires the Nebulosa package for density UMAP plotting.")
}

stopifnot(exists("scRNA1"), inherits(scRNA1, "Seurat"))
dir.create("output", showWarnings = FALSE, recursive = TRUE)

panel_genes <- c("YBX1", "ABCB1", "ABCC1", "ABCC2", "MVP")
transporters <- c("ABCB1", "ABCC1", "ABCC2", "MVP")

find_feature <- function(obj, gene) {
  assays <- names(obj@assays)
  assays <- unique(c(intersect(c("RNA", "SCT"), assays),
                     setdiff(assays, c("RNA", "SCT"))))
  for (assay in assays) {
    rn <- tryCatch(rownames(obj@assays[[assay]]), error = function(e) character(0))
    if (!length(rn)) next
    if (gene %in% rn) return(list(assay = assay, feature = gene))
    hit <- grep(paste0("^", gene, "$"), rn, ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(list(assay = assay, feature = hit[1]))
  }
  stop("Feature not found in any assay: ", gene)
}

feature_info <- setNames(lapply(panel_genes, function(g) find_feature(scRNA1, g)), panel_genes)

## Fig 3F: individual ABC transporter density UMAPs.
density_plots <- lapply(transporters, function(gene) {
  DefaultAssay(scRNA1) <- feature_info[[gene]]$assay
  p <- Nebulosa::plot_density(
    scRNA1,
    features = feature_info[[gene]]$feature,
    reduction = "umap"
  ) + ggtitle(gene)

  ggsave(file.path("output", paste0("FeatureDensity_", gene, "_umap.pdf")),
         p, width = 4, height = 3.6)
  ggsave(file.path("output", paste0("FeatureDensity_", gene, "_umap.png")),
         p, width = 4, height = 3.6, dpi = 300, bg = "white")
  p
})

p_fig3f <- wrap_plots(density_plots, nrow = 1)
ggsave("output/Fig3F_ABC_transporters_density_UMAP.pdf",
       p_fig3f, width = 12, height = 3.6)
ggsave("output/Fig3F_ABC_transporters_density_UMAP.png",
       p_fig3f, width = 12, height = 3.6, dpi = 300, bg = "white")

## Fig S3E-S3H: scaled group means for YB-1 and ABC transporters.
assay_use <- feature_info[["YBX1"]]$assay
DefaultAssay(scRNA1) <- assay_use
expr <- GetAssayData(scRNA1, assay = assay_use, slot = "data")
expr <- as.matrix(expr)
feature_names <- vapply(feature_info, `[[`, character(1), "feature")

make_group_heatmap <- function(group_col, group_levels, outfile, column_labels = group_levels) {
  if (!group_col %in% colnames(scRNA1@meta.data)) {
    stop("Missing metadata column required for heatmap: ", group_col)
  }

  groups <- factor(scRNA1@meta.data[[group_col]], levels = group_levels)
  keep <- !is.na(groups)
  if (!any(keep)) stop("No cells retained for grouping column: ", group_col)

  avg <- sapply(group_levels, function(level) {
    cells <- rownames(scRNA1@meta.data)[keep & groups == level]
    if (!length(cells)) stop("No cells found for ", group_col, " = ", level)
    rowMeans(expr[feature_names, cells, drop = FALSE], na.rm = TRUE)
  })
  colnames(avg) <- column_labels
  rownames(avg) <- c("YB-1", "ABCB1", "ABCC1", "ABCC2", "MVP")

  z <- t(scale(t(avg)))
  z[!is.finite(z)] <- 0

  pheatmap::pheatmap(
    z,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("#2b6cb0", "white", "#b2182b"))(101),
    breaks = seq(-2, 2, length.out = 102),
    border_color = "grey35",
    fontsize = 10,
    filename = file.path("output", outfile),
    width = 2.4,
    height = 3.4
  )
}

make_group_heatmap(
  group_col = "sample_group",
  group_levels = c("Tumor", "Adjacent"),
  outfile = "Heatmap_YBX1_ABC_Tumor_vs_Adjacent.pdf"
)

make_group_heatmap(
  group_col = "group_copykat",
  group_levels = c("aneuploid", "diploid"),
  outfile = "Heatmap_YBX1_ABC_aneuploid_vs_diploid.pdf"
)

make_group_heatmap(
  group_col = "Cisplatin_pred_group",
  group_levels = c("Predicted sensitive", "Predicted resistant"),
  column_labels = c("CIS Sensitive", "CIS Resistant"),
  outfile = "Heatmap_YBX1_ABC_Cisplatin.pdf"
)

make_group_heatmap(
  group_col = "Gemcitabine_pred_group",
  group_levels = c("Predicted sensitive", "Predicted resistant"),
  column_labels = c("GEM Sensitive", "GEM Resistant"),
  outfile = "Heatmap_YBX1_ABC_Gemcitabine.pdf"
)

message("Done: Fig 3F and Fig S3E-S3H ABC transporter panels saved under output/.")
