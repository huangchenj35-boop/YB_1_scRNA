## ============================================================
## Fig 1E — CopyKAT copy-number status UMAP (epithelial subset)
## ------------------------------------------------------------
## Input : scRNA_sub  (loaded by 00_load_data.R) — 12290 epithelial cells on
##         their own re-embedded UMAP; meta col group_copykat. SAME subset
##         embedding as Fig 1C / Fig 1F.
## Output: results/Fig1E_copykat_UMAP.pdf / .png
## Match : paper colours aneuploid=blue, diploid=red; legend on top.
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

if (!exists("RESULTS")) {
  RESULTS <- Sys.getenv("YB1_RESULTS", unset = "../results")
  dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
}

stopifnot(exists("scRNA_sub"), inherits(scRNA_sub, "Seurat"),
          "group_copykat" %in% colnames(scRNA_sub@meta.data))

p <- DimPlot(
  scRNA_sub,
  reduction = "umap",
  group.by  = "group_copykat",
  cols      = c("aneuploid" = "#377EB8", "diploid" = "#E41A1C"),  # paper: aneuploid blue, diploid red
  label     = FALSE,
  pt.size   = 0.3
) + ggtitle(NULL) + theme(legend.position = "top")

ggsave(file.path(RESULTS, "Fig1E_copykat_UMAP.pdf"), plot = p, width = 7, height = 6)
ggsave(file.path(RESULTS, "Fig1E_copykat_UMAP.png"), plot = p, width = 7, height = 6, dpi = 300, bg = "white")
