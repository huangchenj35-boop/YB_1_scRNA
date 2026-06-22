## ============================================================
## Cell Type Annotation and Marker Validation
## Figures: Fig 1A (cell annotation UMAP), Fig S1B (marker violin plots)
## Dataset: GSE138709 (5 tumor + 3 adjacent iCCA samples)
## ============================================================
## Prerequisites:
##   - scRNA1: integrated Seurat object from 00_data_preprocessing_for_FigS1A.R
##   - seurat_clusters should contain 21 clusters for Fig S1A (0-20)
## Output:
##   - scRNA1 with article-standard curate_v1 labels in metadata
##   - output/Fig1A_cell_annotation_UMAP.pdf/.png
##   - output/FigS1B_celltype_marker_violin.pdf/.png
##   - output/celltype_annotation_audit.csv
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
set.seed(1234)

if (!exists("scRNA1") || !inherits(scRNA1, "Seurat")) {
  scRNA1 <- readRDS("output/scRNA1_preprocessed.rds")
}
dir.create("output", showWarnings = FALSE, recursive = TRUE)

## Article-standard labels from PDF Fig 1A / Fig S1B.
expected_celltypes <- c(
  "B cell", "CD8+ T cell", "Cholangiocytes", "Dendritic cells",
  "Endothelial cell", "Fibroblast", "Hepatocytes", "moDCs",
  "Mono/Macro", "NK cells", "Plasma cell",
  "Plasmacytoid dendritic cells", "T cells"
)

## Normalize older labels from the supplied scripts to the exact names in the PDF.
normalize_celltype <- function(x) {
  x <- as.character(x)
  recode <- c(
    "B cells" = "B cell",
    "CD8+ Tex cell" = "CD8+ T cell",
    "CD8+ T cells" = "CD8+ T cell",
    "CD8+ T cell" = "CD8+ T cell",
    "T.cell" = "T cells",
    "T cell" = "T cells",
    "Dendritic cell" = "Dendritic cells",
    "Dendritic cells" = "Dendritic cells",
    "Endothelial cells" = "Endothelial cell",
    "endo" = "Endothelial cell",
    "Fibroblasts" = "Fibroblast",
    "fibr" = "Fibroblast",
    "Macrophage" = "Mono/Macro",
    "Monocyte" = "Mono/Macro",
    "Mono/Macrophage" = "Mono/Macro",
    "NK cell" = "NK cells",
    "Plasma cells" = "Plasma cell",
    "Plasmacytoid dendritic cell" = "Plasmacytoid dendritic cells"
  )
  idx <- x %in% names(recode)
  x[idx] <- recode[x[idx]]
  x
}

## Prefer an existing manually curated annotation if present. If absent, use the
## compact fallback mapping from the supplied V2 annotation script, then audit it.
if ("curate_v1" %in% colnames(scRNA1@meta.data)) {
  scRNA1$curate_v1 <- normalize_celltype(scRNA1$curate_v1)
} else {
  cluster_annotation <- c(
    "0"  = "Hepatocytes",
    "1"  = "T cells",
    "2"  = "Mono/Macro",
    "3"  = "Fibroblast",
    "4"  = "Endothelial cell",
    "5"  = "B cell",
    "6"  = "Cholangiocytes",
    "7"  = "NK cells",
    "8"  = "CD8+ T cell",
    "9"  = "Dendritic cells",
    "10" = "moDCs",
    "11" = "Plasma cell",
    "12" = "Plasmacytoid dendritic cells"
  )
  clusters <- as.character(scRNA1@meta.data$seurat_clusters)
  scRNA1$curate_v1 <- unname(cluster_annotation[clusters])
}

scRNA1$curate_v1 <- factor(normalize_celltype(scRNA1$curate_v1),
                           levels = expected_celltypes)

## Strict audit against the article text and figure legends.
cluster_ids <- sort(unique(as.character(scRNA1@meta.data$seurat_clusters)))
audit <- data.frame(
  item = c(
    "Fig S1A cluster count",
    "Fig 1A cell-type labels present",
    "Fig 1A unexpected cell-type labels",
    "Fig 1A cells with unmapped curate_v1"
  ),
  expected = c(
    "21 clusters (0-20)",
    paste(expected_celltypes, collapse = "; "),
    "none",
    "0"
  ),
  observed = c(
    paste0(length(cluster_ids), " clusters: ", paste(cluster_ids, collapse = ", ")),
    paste(sort(unique(as.character(na.omit(scRNA1$curate_v1)))), collapse = "; "),
    paste(setdiff(sort(unique(as.character(na.omit(scRNA1$curate_v1)))), expected_celltypes),
          collapse = "; "),
    sum(is.na(scRNA1$curate_v1))
  ),
  stringsAsFactors = FALSE
)
write.csv(audit, "output/celltype_annotation_audit.csv", row.names = FALSE)

if (length(cluster_ids) != 21) {
  stop("PDF Fig S1A reports 21 clusters (0-20), but code observed ",
       length(cluster_ids), " clusters: ", paste(cluster_ids, collapse = ", "),
       ". Check preprocessing resolution before generating article figures.")
}
missing_celltypes <- setdiff(expected_celltypes, unique(as.character(na.omit(scRNA1$curate_v1))))
if (length(missing_celltypes)) {
  stop("Missing article cell-type labels in curate_v1: ",
       paste(missing_celltypes, collapse = ", "),
       ". Update cluster-to-celltype annotation before generating Fig 1A/S1B.")
}
if (anyNA(scRNA1$curate_v1)) {
  unmapped <- sort(unique(as.character(scRNA1@meta.data$seurat_clusters[is.na(scRNA1$curate_v1)])))
  stop("Some clusters/cells are not mapped to article cell-type labels: ",
       paste(unmapped, collapse = ", "),
       ". Preserve the original annotation logic, but update its label mapping before generating Fig 1A/S1B.")
}

## Fig 1A: article cell-type annotation UMAP.
p_umap_annot <- DimPlot(
  scRNA1,
  group.by   = "curate_v1",
  label      = TRUE,
  label.size = 3,
  reduction  = "umap"
) + ggtitle("Cell Annotation")

ggsave("output/Fig1A_cell_annotation_UMAP.pdf", p_umap_annot, width = 10, height = 7)
ggsave("output/Fig1A_cell_annotation_UMAP.png", p_umap_annot, width = 10, height = 7,
       dpi = 300, bg = "white")

## Fig S1B: marker violin plot using the markers named in the PDF text/legend.
canonical_markers <- list(
  "B cell"                       = c("CD19", "MS4A1", "CD79A", "CD79B"),
  "CD8+ T cell"                  = c("CD8A", "CD8B", "GZMB", "IFNG"),
  "Cholangiocytes"               = c("SOX9", "KRT19", "EPCAM"),
  "Dendritic cells"              = c("CLEC9A", "XCR1", "BATF3"),
  "Endothelial cell"             = c("VWF", "PECAM1", "ESAM"),
  "Fibroblast"                   = c("COL1A1", "COL6A1", "MYL9"),
  "Hepatocytes"                  = c("ALB", "CYP3A4", "APOA1"),
  "moDCs"                        = c("FCGR2A", "CD86", "CD14"),
  "Mono/Macro"                   = c("FCGR3A", "S100A8", "S100A9"),
  "NK cells"                     = c("CD68", "CSF1R", "IFNG", "KLRC1", "KLRD1", "KLRB1"),
  "Plasma cell"                  = c("IGJ"),
  "Plasmacytoid dendritic cells" = c("SLAMF7", "LILRA4", "IL3RA", "TCF4", "CLEC4C", "IRF7"),
  "T cells"                      = c("CD3D", "CD3E")
)

marker_genes <- unique(unlist(canonical_markers))
marker_genes <- marker_genes[marker_genes %in% rownames(scRNA1)]
if (!length(marker_genes)) {
  stop("None of the article marker genes were found in the Seurat object.")
}

p_vln <- VlnPlot(
  scRNA1,
  features = marker_genes,
  group.by = "curate_v1",
  pt.size  = 0,
  stack    = TRUE,
  flip     = TRUE
) + NoLegend()

ggsave("output/FigS1B_celltype_marker_violin.pdf", p_vln, width = 14, height = 16)
ggsave("output/FigS1B_celltype_marker_violin.png", p_vln, width = 14, height = 16,
       dpi = 300, bg = "white")

saveRDS(scRNA1, file = "output/scRNA1_annotated.rds")
message("Done: Fig 1A/S1B annotation audited and saved.")
