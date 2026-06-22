## ============================================================
## GSVA / ssGSEA: YBX1 Target Gene Activity Scoring
## Figures: Fig 3C top
## Dataset: GSE138709 (5 tumor + 3 adjacent iCCA samples)
## ============================================================
## Prerequisites:
##   - scRNA1: Seurat object
##   - SCENIC regulon target table with columns TF and gene
## Output:
##   - scRNA1 with GSVA_YBX1_targets_ssgsea metadata column
##   - output/FeaturePlot_GSVA_YBX1_targets_ssgsea_<reduction>.pdf/.png
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(GSVA)
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
set.seed(1234)
dir.create("output", showWarnings = FALSE, recursive = TRUE)

if (!exists("scRNA1") || !inherits(scRNA1, "Seurat")) {
  candidates <- c("output/scRNAauc.rds", "output/scRNA1_annotated.rds", "output/scRNA1_preprocessed.rds")
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) {
    stop("scRNA1 was not found in memory and no expected RDS file exists in output/.")
  }
  scRNA1 <- readRDS(hit)
}

tf_root <- function(x) {
  y <- sub("\\s*\\(\\d+g\\)$", "", x)
  y <- sub("_[0-9]+g$", "", y, ignore.case = TRUE)
  y <- sub("_extended$", "", y, ignore.case = TRUE)
  y
}

find_target_file <- function() {
  candidates <- c(
    "output/SCENIC_regulon_target_table_Cisplatin.tsv",
    "output/SCENIC_regulon_target_table.tsv",
    "/home/sodajay/scRNA_TAO/GSE138709/output/SCENIC_regulon_target_table_Cisplatin.tsv"
  )
  hit <- candidates[file.exists(candidates)][1]
  if (!is.na(hit)) return(hit)

  cand <- c(
    list.files("output", pattern = "SCENIC_regulon_target_table.*\\.tsv$",
               full.names = TRUE),
    list.files("int", pattern = "SCENIC_regulon_target_table.*\\.tsv$",
               full.names = TRUE)
  )
  if (length(cand)) return(cand[1])
  stop("No SCENIC regulon target table was found.")
}

target_file <- find_target_file()
targets_df <- data.table::fread(target_file)
if (!all(c("TF", "gene") %in% colnames(targets_df))) {
  stop("Target table must contain columns named TF and gene: ", target_file)
}

ybx1_targets <- targets_df %>%
  mutate(TF_root = tf_root(TF)) %>%
  filter(TF_root == "YBX1") %>%
  pull(gene) %>%
  unique() %>%
  na.omit() %>%
  as.character()

if (!length(ybx1_targets)) stop("No YBX1 target genes were found in target table.")

assay_use <- if ("RNA" %in% names(scRNA1@assays)) "RNA" else DefaultAssay(scRNA1)
DefaultAssay(scRNA1) <- assay_use

expr_data <- tryCatch(GetAssayData(scRNA1, assay = assay_use, slot = "data"),
                      error = function(e) NULL)
use_data <- !is.null(expr_data) && sum(expr_data != 0) > 0
expr_mat <- if (use_data) expr_data else GetAssayData(scRNA1, assay = assay_use, slot = "counts")
kcdf_use <- if (use_data) "Gaussian" else "Poisson"

genes_use <- intersect(rownames(expr_mat), ybx1_targets)
if (length(genes_use) < 3) {
  warning("Fewer than 3 YBX1 target genes overlap the expression matrix: ", length(genes_use))
}
gene_sets <- list(YBX1_targets = genes_use)

run_ssgsea <- function(expr, gene_sets, kcdf) {
  expr <- as.matrix(expr)
  if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    param <- GSVA::ssgseaParam(expr, gene_sets, normalize = TRUE)
    return(GSVA::gsva(param, verbose = FALSE))
  }
  GSVA::gsva(
    expr = expr,
    gset.idx.list = gene_sets,
    method = "ssgsea",
    kcdf = kcdf,
    abs.ranking = TRUE,
    min.sz = 1,
    max.sz = Inf,
    parallel.sz = 1,
    verbose = FALSE
  )
}

scores <- run_ssgsea(expr_mat, gene_sets, kcdf_use)
ybx1_ssgsea <- as.numeric(scores["YBX1_targets", colnames(scRNA1), drop = TRUE])

z_score <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sdv <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) rep(0, length(x)) else (x - mu) / sdv
}

scRNA1$GSVA_YBX1_targets_ssgsea <- ybx1_ssgsea
scRNA1$GSVA_YBX1_targets_ssgsea_z <- z_score(ybx1_ssgsea)

red_use <- if ("umap" %in% names(scRNA1@reductions)) {
  "umap"
} else if ("tsne" %in% names(scRNA1@reductions)) {
  "tsne"
} else {
  stop("scRNA1 has neither UMAP nor tSNE reduction.")
}

p_raw <- FeaturePlot(scRNA1, features = "GSVA_YBX1_targets_ssgsea",
                     label = FALSE, reduction = red_use) +
  ggtitle("YB-1 target ssGSEA")
p_z <- FeaturePlot(scRNA1, features = "GSVA_YBX1_targets_ssgsea_z",
                   label = FALSE, reduction = red_use) +
  ggtitle("YB-1 target ssGSEA z-score")

ggsave(file.path("output", paste0("FeaturePlot_GSVA_YBX1_targets_ssgsea_", red_use, ".pdf")),
       p_raw, width = 7, height = 6)
ggsave(file.path("output", paste0("FeaturePlot_GSVA_YBX1_targets_ssgsea_", red_use, ".png")),
       p_raw, width = 7, height = 6, dpi = 300, bg = "white")
ggsave(file.path("output", paste0("FeaturePlot_GSVA_YBX1_targets_ssgsea_z_", red_use, ".pdf")),
       p_z, width = 7, height = 6)
ggsave(file.path("output", paste0("FeaturePlot_GSVA_YBX1_targets_ssgsea_z_", red_use, ".png")),
       p_z, width = 7, height = 6, dpi = 300, bg = "white")

saveRDS(scRNA1, file = "output/scRNA1_with_YBX1_GSVA.rds")
message("Done: Fig 3C top GSVA/ssGSEA UMAP saved.")
