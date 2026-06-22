## ============================================================
## Cholangiocyte Lineage Subclustering
## Figures: Fig S1D, Fig S1E, Fig S1F
## Dataset: GSE138709 (5 tumor + 3 adjacent iCCA samples)
## ============================================================
## Prerequisites:
##   - scRNA1.subset: cholangiocyte-lineage Seurat object, or
##   - scRNA1 with curate_v1 annotation from 01_Fig1A_FigS1B_cell_annotation.R
## Output:
##   - output/FigS1D_cholangiocyte_subcluster_UMAP.pdf/.png
##   - output/FigS1E_cholangiocyte_subtype_UMAP.pdf/.png
##   - output/FigS1F_cholangiocyte_subtype_marker_violin.pdf/.png
##   - output/cholangiocyte_subtype_annotation_audit.csv
## Notes:
##   - The cluster-to-subtype logic is kept from the supplied script.
##   - Label wording is aligned to the PDF: LPC, Cholangiocyte, Malignant.
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

dir.create("output", showWarnings = FALSE, recursive = TRUE)

if (exists("scRNA1.subset") && inherits(scRNA1.subset, "Seurat")) {
  obj <- scRNA1.subset
} else if (exists("scRNA1") && inherits(scRNA1, "Seurat")) {
  if (!"curate_v1" %in% colnames(scRNA1@meta.data)) {
    stop("scRNA1.subset is absent and scRNA1 has no curate_v1 annotation.")
  }
  chol_cells <- rownames(scRNA1@meta.data)[scRNA1$curate_v1 == "Cholangiocytes"]
  if (!length(chol_cells)) stop("No Cholangiocytes were found in scRNA1$curate_v1.")
  obj <- subset(scRNA1, cells = chol_cells)
} else {
  candidates <- c("output/scRNA1_annotated.rds", "output/scRNA1_preprocessed.rds")
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) {
    stop("No scRNA1.subset/scRNA1 object was found and no expected RDS file exists.")
  }
  scRNA1 <- readRDS(hit)
  if (!"curate_v1" %in% colnames(scRNA1@meta.data)) {
    stop("Loaded scRNA1 has no curate_v1 annotation.")
  }
  chol_cells <- rownames(scRNA1@meta.data)[scRNA1$curate_v1 == "Cholangiocytes"]
  if (!length(chol_cells)) stop("No Cholangiocytes were found in loaded scRNA1$curate_v1.")
  obj <- subset(scRNA1, cells = chol_cells)
}

cluster2type <- c(
  "0"  = "Malignant",
  "1"  = "Malignant",
  "2"  = "Cholangiocyte",
  "3"  = "Malignant",
  "4"  = "Malignant",
  "5"  = "Cholangiocyte",
  "6"  = "Cholangiocyte",
  "7"  = "Cholangiocyte",
  "8"  = "Cholangiocyte",
  "9"  = "Cholangiocyte",
  "10" = "Malignant",
  "11" = "LPC",
  "12" = "Malignant",
  "13" = "Malignant"
)

clusters <- as.character(obj@meta.data$seurat_clusters)
obj@meta.data$curate_v1 <- unname(cluster2type[clusters])
obj@meta.data$curate_v1 <- factor(
  obj@meta.data$curate_v1,
  levels = c("LPC", "Cholangiocyte", "Malignant")
)

audit <- data.frame(
  item = c("mapped clusters", "unmapped clusters", "subtype labels"),
  expected = c(
    "0-13 according to supplied cluster2type mapping",
    "none",
    "LPC; Cholangiocyte; Malignant"
  ),
  observed = c(
    paste(sort(unique(clusters)), collapse = "; "),
    paste(sort(unique(clusters[is.na(obj@meta.data$curate_v1)])), collapse = "; "),
    paste(levels(obj@meta.data$curate_v1), collapse = "; ")
  ),
  stringsAsFactors = FALSE
)
write.csv(audit, "output/cholangiocyte_subtype_annotation_audit.csv", row.names = FALSE)

if (anyNA(obj@meta.data$curate_v1)) {
  stop("Unmapped cholangiocyte subclusters: ",
       paste(sort(unique(clusters[is.na(obj@meta.data$curate_v1)])), collapse = ", "))
}

red_use <- if ("umap" %in% names(obj@reductions)) {
  "umap"
} else if ("tsne" %in% names(obj@reductions)) {
  "tsne"
} else {
  stop("The cholangiocyte object has neither UMAP nor tSNE reduction.")
}

p_s1d <- DimPlot(
  obj,
  reduction = red_use,
  group.by = "seurat_clusters",
  label = TRUE,
  label.size = 4
) + ggtitle("Cholangiocyte subclusters")

p_s1e <- DimPlot(
  obj,
  reduction = red_use,
  group.by = "curate_v1",
  label = TRUE,
  label.size = 4
) + ggtitle("Cholangiocyte lineage annotation")

ggsave(file.path("output", paste0("FigS1D_cholangiocyte_subcluster_", red_use, ".pdf")),
       p_s1d, width = 7, height = 6)
ggsave(file.path("output", paste0("FigS1D_cholangiocyte_subcluster_", red_use, ".png")),
       p_s1d, width = 7, height = 6, dpi = 300, bg = "white")

ggsave(file.path("output", paste0("FigS1E_cholangiocyte_subtype_", red_use, ".pdf")),
       p_s1e, width = 7, height = 6)
ggsave(file.path("output", paste0("FigS1E_cholangiocyte_subtype_", red_use, ".png")),
       p_s1e, width = 7, height = 6, dpi = 300, bg = "white")

DefaultAssay(obj) <- if ("RNA" %in% names(obj@assays)) "RNA" else DefaultAssay(obj)

features <- c(
  "HNF4A", "ALB", "SOX9", "KRT19", "MYC",
  "KI67", "TOP2A", "CENPF", "BIRC5", "UBE2C",
  "CDK1", "CCNB1", "CCNB2", "VIM", "S100A10",
  "S100A11", "TAGLN", "FN1", "SPARC", "COL1A1",
  "COL1A2", "ITGA6", "ITGB1", "ZEB1", "TWIST1"
)
features <- features[features %in% rownames(obj[[DefaultAssay(obj)]])]
if (!length(features)) {
  stop("None of the Fig S1F marker genes were found in the active assay.")
}

vln_df <- FetchData(obj, vars = c(features, "curate_v1")) %>%
  tibble::rownames_to_column("cell") %>%
  tidyr::pivot_longer(
    cols = all_of(features),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  dplyr::mutate(
    gene = factor(gene, levels = features),
    curate_v1 = factor(curate_v1, levels = c("LPC", "Cholangiocyte", "Malignant"))
  )

p_s1f <- ggplot(vln_df, aes(x = curate_v1, y = expression, fill = curate_v1)) +
  geom_violin(scale = "width", linewidth = 0.2, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.2, fill = "white") +
  facet_wrap(~ gene, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c(
    "LPC" = "#6A7FDB",
    "Cholangiocyte" = "#52A97F",
    "Malignant" = "#D65A4A"
  )) +
  labs(x = NULL, y = "Expression") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("output/FigS1F_cholangiocyte_subtype_marker_violin.pdf",
       p_s1f, width = 13, height = 10)
ggsave("output/FigS1F_cholangiocyte_subtype_marker_violin.png",
       p_s1f, width = 13, height = 10, dpi = 300, bg = "white")

saveRDS(obj, file = "output/scRNA1_cholangiocyte_subtypes.rds")
message("Done: Fig S1D/S1E/S1F cholangiocyte panels saved.")
