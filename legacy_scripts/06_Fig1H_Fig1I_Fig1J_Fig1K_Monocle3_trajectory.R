## ============================================================
## Pseudotime Trajectory Analysis (Monocle3)
## Figures: Fig 1H, Fig 1I, Fig 1J, Fig 1K
## Dataset: GSE138709 (5 tumor + 3 adjacent iCCA samples)
## ============================================================
## Prerequisites:
##   - scRNA1.subset trajectory input, or scRNA1
##   - group_copykat metadata from 04_Fig1E_FigS2_inferCNV_CopyKAT.R
## Output:
##   - output/monocle3_CD8T_pseudotime.pdf/.png          (Fig 1H)
##   - output/monocle3_copykat.pdf/.png                  (Fig 1I)
##   - output/pseudotime_trajectory_by_state.pdf/.png    (Fig 1J)
##   - output/YBX1_in_pseudotime.pdf/.png                (Fig 1K)
##   - output/YBX1_pseudotime_smooth.pdf/.png            (Fig 1K)
## ============================================================

suppressPackageStartupMessages({
  library(monocle3)
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
})

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
dir.create("output", showWarnings = FALSE, recursive = TRUE)

if (exists("scRNA1.subset") && inherits(scRNA1.subset, "Seurat")) {
  obj <- scRNA1.subset
} else if (exists("scRNA1") && inherits(scRNA1, "Seurat")) {
  obj <- scRNA1
} else {
  candidates <- c("output/scRNA1_annotated.rds", "output/scRNA1_preprocessed.rds")
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) {
    stop("No scRNA1.subset/scRNA1 object was found and no expected RDS exists.")
  }
  obj <- readRDS(hit)
}

if (!"group_copykat" %in% colnames(obj@meta.data)) {
  stop("group_copykat metadata is missing. Run 04_Fig1E_FigS2_inferCNV_CopyKAT.R first.")
}

DefaultAssay(obj) <- if ("RNA" %in% names(obj@assays)) "RNA" else DefaultAssay(obj)
counts <- GetAssayData(obj, assay = DefaultAssay(obj), slot = "counts")
counts <- counts[!duplicated(rownames(counts)), , drop = FALSE]
counts <- as(counts, "CsparseMatrix")

cell_metadata <- obj@meta.data
colnames(cell_metadata) <- make.unique(colnames(cell_metadata))
gene_annotation <- data.frame(
  gene_short_name = rownames(counts),
  row.names = rownames(counts),
  stringsAsFactors = FALSE
)

cds <- new_cell_data_set(
  counts,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

cds <- preprocess_cds(cds, num_dim = min(100, ncol(cds) - 1))
if ("orig.ident" %in% colnames(colData(cds))) {
  cds <- align_cds(cds, num_dim = min(100, ncol(cds) - 1), alignment_group = "orig.ident")
}
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds, resolution = 1e-5, reduction_method = "UMAP")
cds <- learn_graph(cds)

get_earliest_principal_node <- function(cds, time_bin = "diploid") {
  cell_ids <- which(colData(cds)[, "group_copykat"] == time_bin)
  if (!length(cell_ids)) return(NULL)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), , drop = FALSE])
  igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
  ]
}

root_node <- get_earliest_principal_node(cds, "diploid")
if (!is.null(root_node)) {
  cds <- order_cells(cds, root_pr_nodes = root_node)
} else {
  warning("No diploid root node was found; ordering cells without explicit root.")
  cds <- order_cells(cds)
}

p_h <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 1.5,
  group_label_size = 4,
  cell_size = 1.2
)

p_i <- plot_cells(
  cds,
  color_cells_by = "group_copykat",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  cell_size = 1.2
)

p_j <- plot_cells(
  cds,
  color_cells_by = "cluster",
  label_cell_groups = TRUE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 1.5,
  group_label_size = 4,
  cell_size = 1.2
)

ggsave("output/monocle3_CD8T_pseudotime.pdf", p_h, width = 10, height = 6)
ggsave("output/monocle3_CD8T_pseudotime.png", p_h, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("output/monocle3_copykat.pdf", p_i, width = 10, height = 6)
ggsave("output/monocle3_copykat.png", p_i, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("output/pseudotime_trajectory_by_state.pdf", p_j, width = 10, height = 6)
ggsave("output/pseudotime_trajectory_by_state.png", p_j, width = 10, height = 6, dpi = 300, bg = "white")

if (!"YBX1" %in% rownames(cds)) {
  stop("YBX1 was not found in the Monocle3 CellDataSet.")
}

p_ybx1_pt <- plot_genes_in_pseudotime(
  cds["YBX1", ],
  min_expr = 0.1
) + ggtitle("YBX1 along pseudotime")

df_ybx1 <- tibble(
  pseudotime = pseudotime(cds),
  group_copykat = as.character(colData(cds)$group_copykat),
  YBX1 = as.numeric(monocle3::exprs(cds)["YBX1", ])
) %>%
  tidyr::drop_na(pseudotime)

p_ybx1_smooth <- ggplot(df_ybx1, aes(pseudotime, YBX1, color = group_copykat)) +
  geom_point(size = 0.2, alpha = 0.35) +
  geom_smooth(se = FALSE, linewidth = 0.8) +
  theme_classic() +
  labs(x = "Pseudotime", y = "YBX1 expression", color = NULL)

ggsave("output/YBX1_in_pseudotime.pdf", p_ybx1_pt, width = 8, height = 5)
ggsave("output/YBX1_in_pseudotime.png", p_ybx1_pt, width = 8, height = 5, dpi = 300, bg = "white")
ggsave("output/YBX1_pseudotime_smooth.pdf", p_ybx1_smooth, width = 8, height = 5)
ggsave("output/YBX1_pseudotime_smooth.png", p_ybx1_smooth, width = 8, height = 5, dpi = 300, bg = "white")

saveRDS(cds, "output/monocle3_cds_Fig1H_to_Fig1K.rds")
message("Done: Fig 1H-K Monocle3 trajectory panels saved.")
