## ============================================================
## Fig 3C (top) — GSVA YB-1 target score UMAP (epithelial subset)
## Source : 10_Script42_Fig1C_1F_3B_S3D.R — FeaturePlot of GSVA_YBX1_targets_ssgsea
## Input  : scRNAauc (meta col GSVA_YBX1_targets_ssgsea; reduction umap)
## Output : results/Fig3C_GSVA_YBX1_targets_UMAP.pdf / .png
## ============================================================
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
if (!exists("RESULTS")) { RESULTS <- Sys.getenv("YB1_RESULTS", unset="../results"); dir.create(RESULTS, recursive=TRUE, showWarnings=FALSE) }
stopifnot(exists("scRNAauc"), "GSVA_YBX1_targets_ssgsea" %in% colnames(scRNAauc@meta.data))
p <- FeaturePlot(scRNAauc, features = "GSVA_YBX1_targets_ssgsea", reduction = "umap", order = TRUE) + ggtitle("GSVA YB-1 targets")
ggsave(file.path(RESULTS, "Fig3C_GSVA_YBX1_targets_UMAP.pdf"), p, width = 7, height = 6)
ggsave(file.path(RESULTS, "Fig3C_GSVA_YBX1_targets_UMAP.png"), p, width = 7, height = 6, dpi = 300, bg = "white")
