## ============================================================
## Fig 3B (top) — YB-1 regulon AUC UMAP (epithelial subset)
## Source : 10_Script42_Fig1C_1F_3B_S3D.R — FeaturePlot of YBX1_extended_907g
## Input  : scRNAauc (meta col YBX1_extended_907g; reduction umap)
## Output : results/Fig3B_YBX1_regulon_AUC_UMAP.pdf / .png
## ============================================================
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
if (!exists("RESULTS")) { RESULTS <- Sys.getenv("YB1_RESULTS", unset="../results"); dir.create(RESULTS, recursive=TRUE, showWarnings=FALSE) }
stopifnot(exists("scRNAauc"), "YBX1_extended_907g" %in% colnames(scRNAauc@meta.data))
p <- FeaturePlot(scRNAauc, features = "YBX1_extended_907g", reduction = "umap", order = TRUE) + ggtitle("YB-1 regulonAUC")
ggsave(file.path(RESULTS, "Fig3B_YBX1_regulon_AUC_UMAP.pdf"), p, width = 7, height = 6)
ggsave(file.path(RESULTS, "Fig3B_YBX1_regulon_AUC_UMAP.png"), p, width = 7, height = 6, dpi = 300, bg = "white")
