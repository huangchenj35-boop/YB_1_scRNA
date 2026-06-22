## ============================================================
## inferCNV + CopyKAT: Malignant Cell Identification
## Figures: Fig S2 (inferCNV heatmap), Fig 1E (CopyKAT UMAP)
## Dataset: GSE138709 (5 tumor + 3 adjacent iCCA samples)
## ============================================================
## Prerequisites:
##   - scRNA1: Seurat object with cell type annotation (curate_v1)
##             produced by 01_Fig1A_FigS1B_cell_annotation.R
##   - GSE138709/cellAnnotations.txt  (barcode → sample mapping)
##   - GSE138709/geneFile.txt         (gene → chr coordinates, hg38)
## Output:
##   - infercnv.png / infercnv.pdf    (Fig S2 heatmap)
##   - scRNA1 with group_copykat metadata (diploid/aneuploid)
##   - group_copykat_UMAP.pdf          (Fig 1E)
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(infercnv)
  library(copykat)
  library(ggplot2)
})

# ------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------
data_dir   <- "GSE138709"
annot_file <- file.path(data_dir, "cellAnnotations.txt")
gene_file  <- file.path(data_dir, "geneFile.txt")
outdir     <- "output"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ------------------------------------------------------------------
# 1. Prepare inferCNV inputs from scRNA1
# ------------------------------------------------------------------
# Raw count matrix (genes x cells) – inferCNV requires integer counts
raw_counts <- GetAssayData(scRNA1, assay = "RNA", slot = "counts")

# Cell annotation: T cells from adjacent tissue = diploid reference
# cellAnnotations.txt format: barcode \t cell_type
# We use the curate_v1 labels where "T cells" in adjacent samples are reference
cell_annot <- data.frame(
  barcode   = colnames(scRNA1),
  cell_type = scRNA1$curate_v1,
  stringsAsFactors = FALSE
)
# Write annotation file (required by infercnv::CreateInfercnvObject)
annot_out <- file.path(data_dir, "cellAnnotations_infercnv.txt")
write.table(cell_annot, annot_out,
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# ------------------------------------------------------------------
# 2. Run inferCNV
# ------------------------------------------------------------------
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = raw_counts,
  annotations_file  = annot_out,
  delim             = "\t",
  gene_order_file   = gene_file,
  ref_group_names   = c("T cells")   # diploid reference: T cells from adjacent
)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff                    = 0.1,    # minimum mean expression for a gene
  out_dir                   = data_dir,
  cluster_by_groups         = TRUE,
  denoise                   = TRUE,   # default settings including denoising
  HMM                       = TRUE,
  num_threads               = 4,
  output_format             = "png"
)

message("✅ inferCNV finished. Heatmap saved to: ", data_dir)

# ------------------------------------------------------------------
# 3. Classify cells as diploid / aneuploid using CopyKAT
# ------------------------------------------------------------------
# CopyKAT uses raw count matrix of cholangiocyte-lineage cells
# (malignant + cholangiocyte + LPC from scRNA1)
chol_cells <- colnames(scRNA1)[scRNA1$curate_v1 %in%
                c("Malignant", "Cholangiocytes", "LPC",
                  "malignant", "cholangiocyte", "Cholangiocyte")]

if (length(chol_cells) == 0) {
  # Fallback: use all non-T, non-B immune cells
  chol_cells <- colnames(scRNA1)[!scRNA1$curate_v1 %in%
                  c("T cells", "CD8+ T cell", "B cell",
                    "NK cells", "Plasma cell",
                    "Dendritic cells", "Plasmacytoid dendritic cells",
                    "moDCs", "Mono/Macro")]
}

chol_counts <- as.matrix(
  GetAssayData(scRNA1[, chol_cells], assay = "RNA", slot = "counts")
)

# Run CopyKAT (normal cells = adjacent T cells as reference)
t_adj_cells <- colnames(scRNA1)[
  scRNA1$curate_v1 %in% c("T cells", "CD8+ T cell") &
  scRNA1$sample_type == "Adjacent"
]
if (length(t_adj_cells) == 0) {
  t_adj_cells <- colnames(scRNA1)[scRNA1$curate_v1 %in% c("T cells", "CD8+ T cell")]
}

copykat_result <- copykat(
  rawmat        = chol_counts,
  id.type       = "S",            # Seurat-style barcodes
  cell.line     = "no",
  distance      = "euclidean",
  norm.cell.names = t_adj_cells,  # diploid normal reference
  output.seg    = FALSE,
  plot.genes    = TRUE,
  genome        = "hg20",
  n.cores       = 4
)

# ------------------------------------------------------------------
# 4. Add group_copykat to scRNA1 metadata
# ------------------------------------------------------------------
pred <- copykat_result$prediction
# copykat labels: "aneuploid", "diploid", "not.defined"
pred_df <- data.frame(
  barcode       = pred$cell.names,
  group_copykat = pred$copykat.pred,
  stringsAsFactors = FALSE
)
pred_df$group_copykat[pred_df$group_copykat == "not.defined"] <- NA

# Map to all cells in scRNA1
meta_add <- pred_df$group_copykat
names(meta_add) <- pred_df$barcode
scRNA1 <- AddMetaData(scRNA1, metadata = meta_add, col.name = "group_copykat")

message("✅ group_copykat added to scRNA1 metadata.")
message("   Table: ")
print(table(scRNA1$group_copykat, useNA = "ifany"))

# ------------------------------------------------------------------
# 5. UMAP visualization colored by group_copykat (Fig 1E)
# ------------------------------------------------------------------
p_umap <- DimPlot(
  scRNA1,
  reduction  = "umap",
  group.by   = "group_copykat",
  cols       = c("aneuploid" = "#E41A1C", "diploid" = "#377EB8"),
  label      = FALSE,
  pt.size    = 0.3
) +
  ggtitle("Copy Number Status (CopyKAT)") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(outdir, "group_copykat_UMAP.pdf"),
       plot = p_umap, width = 7, height = 6)
ggsave(file.path(outdir, "group_copykat_UMAP.png"),
       plot = p_umap, width = 7, height = 6, dpi = 300, bg = "white")

message("✅ Fig 1E saved: group_copykat_UMAP.pdf")

# ------------------------------------------------------------------
# 6. Save updated Seurat object
# ------------------------------------------------------------------
saveRDS(scRNA1, "scRNA1_with_copykat.rds")
message("✅ scRNA1 with group_copykat saved to scRNA1_with_copykat.rds")
