## ============================================================
## Fig S1E — Cholangiocyte subtype UMAP (LPC / Cholangiocyte / Malignant)
## ------------------------------------------------------------
## Source : 03_Script82_cholangiocyte_subclustering.R lines 13-45
##          (cluster -> subtype mapping cluster2type), plotted with the
##          DimPlot pattern from 02_Script16 line 405.
## Input  : scRNA_sub (loaded by 00_load_data.R; meta col seurat_clusters)
## Output : results/FigS1E_subtype_UMAP.pdf / .png
##          Subtype labels: LPC (= Hepatic progenitor cell), Cholangiocyte,
##          Malignant (= Malignant cholangiocyte).
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

## ---- cluster -> subtype mapping (verbatim from 03_Script82) ----
cluster2type <- c(
  `0`  = "Malignant cholangiocyte",
  `1`  = "Malignant cholangiocyte",
  `2`  = "Cholangiocyte",
  `3`  = "Malignant cholangiocyte",
  `4`  = "Malignant cholangiocyte",
  `5`  = "Cholangiocyte",
  `6`  = "Cholangiocyte",
  `7`  = "Cholangiocyte",
  `8`  = "Cholangiocyte",
  `9`  = "Cholangiocyte",
  `10` = "Malignant cholangiocyte",
  `11` = "Hepatic progenitor cell",
  `12` = "Malignant cholangiocyte",
  `13` = "Malignant cholangiocyte"
)
clusters <- as.character(scRNA_sub@meta.data$seurat_clusters)
subtype  <- unname(cluster2type[clusters])

## display labels to match the paper legend
subtype[subtype == "Hepatic progenitor cell"] <- "LPC"
subtype[subtype == "Malignant cholangiocyte"] <- "Malignant"
## 与文章 S1E 一致：默认 hue 配色，图例/因子顺序 Cholangiocyte / LPC / Malignant
## （3 色默认 hue 字母序 -> Cholangiocyte=salmon, LPC=green, Malignant=blue）
scRNA_sub$subtype <- factor(subtype, levels = c("Cholangiocyte", "LPC", "Malignant"))

p <- DimPlot(scRNA_sub, group.by = "subtype", label = FALSE, reduction = "umap")

ggsave(file.path(RESULTS, "FigS1E_subtype_UMAP.pdf"), plot = p, width = 7, height = 6)
ggsave(file.path(RESULTS, "FigS1E_subtype_UMAP.png"), plot = p, width = 7, height = 6, dpi = 300, bg = "white")
