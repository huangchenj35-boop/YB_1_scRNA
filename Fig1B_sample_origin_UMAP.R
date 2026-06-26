## ============================================================
## Fig 1B — Sample-origin UMAP (Adjacent vs Tumor)
## ------------------------------------------------------------
## Input : scRNA1   (loaded by 00_load_data.R; meta col `sample_group`)
## Output: results/Fig1B_sample_type_UMAP.pdf / .png
## Note  : plotting parameters are unchanged from the original script;
##         only headers/library calls were added and dead code removed.
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# allow standalone use if 00_load_data.R was not sourced
if (!exists("RESULTS")) {
  RESULTS <- Sys.getenv("YB1_RESULTS", unset = "../results")
  dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
}

# 与文章原图一致：Seurat 默认双色 (paper Fig1B: Adjacent=salmon, Tumor=teal)
pal_contrast <- c("Adjacent" = "#F8766D",   # 肉粉 salmon
                  "Tumor"    = "#00BFC4")   # 青 teal

stopifnot(exists("scRNA1"), inherits(scRNA1, "Seurat"))
scRNA1$sample_type <- factor(scRNA1$sample_group, levels = c("Adjacent", "Tumor"))

p <- DimPlot(
  scRNA1,
  group.by   = "sample_type",
  label      = FALSE,
  label.size = 5,
  reduction  = "umap",
  cols       = unname(pal_contrast[levels(scRNA1$sample_type)])
) + ggtitle("Sample group")
ggsave(file.path(RESULTS, "Fig1B_sample_type_UMAP.pdf"), plot = p, width = 5, height = 4)
ggsave(file.path(RESULTS, "Fig1B_sample_type_UMAP.png"), plot = p, width = 5, height = 4, dpi = 300, bg = "white")
