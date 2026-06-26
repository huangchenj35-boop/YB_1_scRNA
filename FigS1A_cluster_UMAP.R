## ============================================================
## Fig S1A — Seurat cluster UMAP (full object)
## ------------------------------------------------------------
## Source : 01_Script3_preprocessing.R / 02_Script16_cell_annotation.R
##          DimPlot(scRNA1, reduction = "umap", label = TRUE)
## Input  : scRNA1 (loaded by 00_load_data.R; meta col seurat_clusters)
## Output : results/FigS1A_cluster_UMAP.pdf / .png
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

if (!exists("RESULTS")) {
  RESULTS <- Sys.getenv("YB1_RESULTS", unset = "../results")
  dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
}

stopifnot(exists("scRNA1"), "seurat_clusters" %in% colnames(scRNA1@meta.data))

Idents(scRNA1) <- "seurat_clusters"
p <- DimPlot(scRNA1, reduction = "umap", label = TRUE)

ggsave(file.path(RESULTS, "FigS1A_cluster_UMAP.pdf"), plot = p, width = 7, height = 6)
ggsave(file.path(RESULTS, "FigS1A_cluster_UMAP.png"), plot = p, width = 7, height = 6, dpi = 300, bg = "white")
