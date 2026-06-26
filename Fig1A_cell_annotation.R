## ============================================================
## Fig 1A — Cell annotation UMAP (13 cell types)
## ------------------------------------------------------------
## Input : scRNA1   (loaded by 00_load_data.R; meta col curate_v1)
## Output: results/Fig1A_cell_annotation_UMAP.pdf / .png
## Match : paper uses default ggplot/Seurat palette, NO on-plot labels,
##         legend label "CD8+ T cell" (object stores "CD8+ Tex cell"),
##         title "Cell Annotation".
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

if (!exists("RESULTS")) {
  RESULTS <- Sys.getenv("YB1_RESULTS", unset = "../results")
  dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
}

stopifnot(exists("scRNA1"), inherits(scRNA1, "Seurat"))

# display labels to match the paper legend (object level "CD8+ Tex cell" -> "CD8+ T cell")
cv <- as.character(scRNA1$curate_v1)
cv[cv == "CD8+ Tex cell"] <- "CD8+ T cell"
scRNA1$cell_annotation <- factor(cv, levels = sort(unique(cv)))

p <- DimPlot(
  scRNA1,
  group.by  = "cell_annotation",
  label     = FALSE,
  reduction = "umap"
) + ggtitle("Cell Annotation")

ggsave(file.path(RESULTS, "Fig1A_cell_annotation_UMAP.pdf"), plot = p, width = 7, height = 5)
ggsave(file.path(RESULTS, "Fig1A_cell_annotation_UMAP.png"), plot = p, width = 7, height = 5, dpi = 300, bg = "white")
