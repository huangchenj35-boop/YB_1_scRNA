## ============================================================
## Fig S1D — Cholangiocyte (epithelial) sub-cluster UMAP
## ------------------------------------------------------------
## Source : 02_Script16_cell_annotation.R line 408
##          DimPlot(scRNA1.subset, group.by="seurat_clusters",
##                  label=T, label.size=5, reduction='umap')
## Input  : scRNA_sub (loaded by 00_load_data.R) = the epithelial subset
##          (12290 cells) with its own UMAP and seurat_clusters.
## Output : results/FigS1D_subcluster_UMAP.pdf / .png
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

if (!exists("RESULTS")) {
  RESULTS <- Sys.getenv("YB1_RESULTS", unset = "../results")
  dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
}

stopifnot(exists("scRNA_sub"), "seurat_clusters" %in% colnames(scRNA_sub@meta.data))

p <- DimPlot(scRNA_sub, group.by = "seurat_clusters",
             label = TRUE, label.size = 5, reduction = "umap")

ggsave(file.path(RESULTS, "FigS1D_subcluster_UMAP.pdf"), plot = p, width = 7, height = 6)
ggsave(file.path(RESULTS, "FigS1D_subcluster_UMAP.png"), plot = p, width = 7, height = 6, dpi = 300, bg = "white")
