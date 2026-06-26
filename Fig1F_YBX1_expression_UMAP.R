## ============================================================
## Fig 1F — YB-1 expression UMAP (epithelial subset)
## ------------------------------------------------------------
## Input : scRNA_sub (loaded by 00_load_data.R) — epithelial subset, own UMAP;
##         gene YBX1. SAME subset embedding as Fig 1C / Fig 1E.
## Output: results/Fig1F_YBX1_expression_UMAP.pdf / .png
## Match : paper uses a purple gradient (light lavender -> dark indigo), title "YB-1".
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

if (!exists("RESULTS")) {
  RESULTS <- Sys.getenv("YB1_RESULTS", unset = "../results")
  dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
}

stopifnot(exists("scRNA_sub"), inherits(scRNA_sub, "Seurat"))
if (!("YBX1" %in% rownames(scRNA_sub))) DefaultAssay(scRNA_sub) <- "RNA"

p <- FeaturePlot(scRNA_sub, features = "YBX1", reduction = "umap", order = TRUE) +
  scale_color_gradient(low = "#ECE7F4", high = "#3F2C8E", name = NULL) +
  ggtitle("YB-1")

ggsave(file.path(RESULTS, "Fig1F_YBX1_expression_UMAP.pdf"), plot = p, width = 7, height = 6)
ggsave(file.path(RESULTS, "Fig1F_YBX1_expression_UMAP.png"), plot = p, width = 7, height = 6, dpi = 300, bg = "white")
