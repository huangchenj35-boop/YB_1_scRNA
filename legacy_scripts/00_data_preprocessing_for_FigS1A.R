## ============================================================
## scRNA-seq Preprocessing and Integration
## Figures: Fig S1A
## Dataset: GSE138709 (5 tumor + 3 adjacent iCCA samples)
## ============================================================
## Prerequisites:
##   - GEO processed UMI matrices: *_UMI.csv.gz, or
##   - GSE138709_RAW.tar containing *_UMI.csv.gz files, or
##   - 10x CellRanger filtered_feature_bc_matrix folders, or
##   - an existing scRNA1 object / output/scRNA1_preprocessed.rds
## Output:
##   - output/scRNA1_preprocessed.rds
##   - output/FigS1A_cluster_UMAP.pdf/.png
## Notes:
##   - QC thresholds follow the supplied preprocessing logic:
##     nFeature_RNA > 200, nFeature_RNA < 8000,
##     percent.mt < 20, percent.HB < 5, nCount_RNA < 100000.
##   - cluster_resolution defaults to 0.4, matching the supplied integrated
##     cluster workflow. Adjust only if the source object was generated with a
##     different documented resolution.
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
  library(Matrix)
})

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
set.seed(1234)
dir.create("output", showWarnings = FALSE, recursive = TRUE)

if (!exists("cluster_resolution")) cluster_resolution <- 0.4
if (!exists("pc_dims")) pc_dims <- 1:20

infer_sample_group <- function(sample_id) {
  ifelse(
    grepl("tumou?r|tumor|_T$|T$", sample_id, ignore.case = TRUE),
    "Tumor",
    ifelse(grepl("adjacent|normal|_L$|L$", sample_id, ignore.case = TRUE),
           "Adjacent", NA_character_)
  )
}

clean_sample_id <- function(file_path) {
  sample_id <- basename(file_path)
  sample_id <- sub("\\.csv\\.gz$", "", sample_id, ignore.case = TRUE)
  sample_id <- sub("\\.csv$", "", sample_id, ignore.case = TRUE)
  sample_id <- sub("_UMI$", "", sample_id, ignore.case = TRUE)
  sample_id <- sub("^GSM[0-9]+_", "", sample_id)
  sample_id
}

extract_geo_raw_tar <- function() {
  tar_candidates <- c(
    "GSE138709_RAW.tar",
    file.path("GSE138709", "GSE138709_RAW.tar"),
    file.path("GSE138709_RAW", "GSE138709_RAW.tar"),
    file.path("/data", "GSE138709_RAW.tar"),
    file.path("../data", "GSE138709_RAW.tar"),
    file.path("data", "GSE138709_RAW.tar")
  )
  tar_candidates <- tar_candidates[file.exists(tar_candidates)]

  if (length(tar_candidates) == 0) {
    return(invisible(FALSE))
  }

  if (!dir.exists("GSE138709_RAW")) {
    dir.create("GSE138709_RAW", showWarnings = FALSE, recursive = TRUE)
  }

  csv_existing <- list.files(
    "GSE138709_RAW",
    pattern = "_UMI\\.csv(\\.gz)?$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(csv_existing) > 0) {
    return(invisible(TRUE))
  }

  message("Extracting ", tar_candidates[[1]], " to GSE138709_RAW/")
  utils::untar(tar_candidates[[1]], exdir = "GSE138709_RAW")
  invisible(TRUE)
}

discover_umi_csv_files <- function() {
  if (exists("umi_csv_files")) return(umi_csv_files)

  extract_geo_raw_tar()

  roots <- c(
    "GSE138709_RAW",
    "GSE138709",
    "Rawcount",
    "rawcount",
    "data",
    "/data",
    "."
  )
  roots <- unique(roots[dir.exists(roots)])

  hits <- unique(unlist(lapply(roots, function(root) {
    list.files(
      root,
      pattern = "_UMI\\.csv(\\.gz)?$",
      recursive = TRUE,
      full.names = TRUE
    )
  })))

  hits <- hits[!grepl("/output/", hits)]
  hits[order(basename(hits))]
}

discover_10x_dirs <- function() {
  if (exists("raw_dirs")) return(raw_dirs)
  roots <- c("GSE138709", "Rawcount/filtered_feature_bc_matrix", "filtered_feature_bc_matrix", ".")
  roots <- roots[dir.exists(roots)]
  hits <- unique(unlist(lapply(roots, function(root) {
    dirs <- list.dirs(root, recursive = TRUE, full.names = TRUE)
    dirs[file.exists(file.path(dirs, "matrix.mtx")) |
           file.exists(file.path(dirs, "matrix.mtx.gz"))]
  })))
  hits
}

read_umi_csv_counts <- function(file_path) {
  sample_id <- clean_sample_id(file_path)
  message("Reading UMI CSV matrix: ", basename(file_path), " as ", sample_id)

  con <- if (grepl("\\.gz$", file_path, ignore.case = TRUE)) {
    gzfile(file_path, open = "rt")
  } else {
    file(file_path, open = "rt")
  }
  on.exit(close(con), add = TRUE)

  counts_df <- read.csv(
    con,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  counts_mat <- as.matrix(counts_df)
  storage.mode(counts_mat) <- "numeric"
  rownames(counts_mat) <- make.unique(rownames(counts_mat))
  colnames(counts_mat) <- paste(sample_id, colnames(counts_mat), sep = "_")

  Matrix::Matrix(counts_mat, sparse = TRUE)
}

add_qc_metrics <- function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  hb_genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1",
                "HBG2", "HBM", "HBQ1", "HBZ")
  hb_genes <- intersect(hb_genes, rownames(obj))
  obj[["percent.HB"]] <- if (length(hb_genes)) {
    PercentageFeatureSet(obj, features = hb_genes)
  } else {
    rep(0, ncol(obj))
  }
  obj
}

if (exists("scRNA1") && inherits(scRNA1, "Seurat")) {
  message("Using scRNA1 from memory.")
} else if (file.exists("output/scRNA1_preprocessed.rds")) {
  scRNA1 <- readRDS("output/scRNA1_preprocessed.rds")
  message("Loaded existing output/scRNA1_preprocessed.rds.")
} else if (file.exists("output/scRNA1_annotated.rds")) {
  scRNA1 <- readRDS("output/scRNA1_annotated.rds")
  message("Loaded existing output/scRNA1_annotated.rds.")
} else {
  umi_files <- discover_umi_csv_files()
  data_dirs <- discover_10x_dirs()

  seurat_list <- list()

  if (length(umi_files) > 0) {
    sample_ids <- make.unique(vapply(umi_files, clean_sample_id, FUN.VALUE = character(1)))
    seurat_list <- vector("list", length(umi_files))
    names(seurat_list) <- sample_ids

    for (i in seq_along(umi_files)) {
      counts <- read_umi_csv_counts(umi_files[i])
      obj <- CreateSeuratObject(
        counts,
        project = sample_ids[i],
        min.cells = 3,
        min.features = 200
      )
      obj$orig.ident <- sample_ids[i]
      obj$sample_group <- infer_sample_group(sample_ids[i])
      obj <- add_qc_metrics(obj)
      obj <- subset(
        obj,
        subset = nFeature_RNA > 200 &
          nFeature_RNA < 8000 &
          percent.mt < 20 &
          percent.HB < 5 &
          nCount_RNA < 100000
      )
      seurat_list[[i]] <- obj
    }
  } else if (length(data_dirs) > 0) {
    sample_ids <- basename(dirname(data_dirs))
    same_name <- duplicated(sample_ids) | duplicated(sample_ids, fromLast = TRUE)
    sample_ids[same_name] <- basename(data_dirs[same_name])
    sample_ids <- make.unique(sample_ids)

    seurat_list <- vector("list", length(data_dirs))
    names(seurat_list) <- sample_ids

    for (i in seq_along(data_dirs)) {
      counts <- Read10X(data.dir = data_dirs[i])
      obj <- CreateSeuratObject(counts, project = sample_ids[i],
                                min.cells = 3, min.features = 200)
      obj$orig.ident <- sample_ids[i]
      obj$sample_group <- infer_sample_group(sample_ids[i])
      obj <- add_qc_metrics(obj)
      obj <- subset(
        obj,
        subset = nFeature_RNA > 200 &
          nFeature_RNA < 8000 &
          percent.mt < 20 &
          percent.HB < 5 &
          nCount_RNA < 100000
      )
      seurat_list[[i]] <- obj
    }
  } else {
    stop(
      "No input count matrices found. Provide GSE138709_RAW.tar, *_UMI.csv.gz files, ",
      "10x filtered_feature_bc_matrix folders, or output/scRNA1_preprocessed.rds."
    )
  }

  if (length(seurat_list) > 1) {
    seurat_list <- lapply(seurat_list, function(obj) {
      obj <- NormalizeData(obj, verbose = FALSE)
      obj <- FindVariableFeatures(obj, selection.method = "vst",
                                  nfeatures = 2000, verbose = FALSE)
      obj
    })
    anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = pc_dims)
    scRNA1 <- IntegrateData(anchorset = anchors, dims = pc_dims)
    DefaultAssay(scRNA1) <- "integrated"
  } else {
    scRNA1 <- seurat_list[[1]]
    scRNA1 <- NormalizeData(scRNA1, verbose = FALSE)
    scRNA1 <- FindVariableFeatures(scRNA1, selection.method = "vst",
                                   nfeatures = 2000, verbose = FALSE)
    DefaultAssay(scRNA1) <- "RNA"
  }

  scRNA1 <- ScaleData(scRNA1, verbose = FALSE)
  scRNA1 <- RunPCA(scRNA1, npcs = max(pc_dims), verbose = FALSE)
  scRNA1 <- FindNeighbors(scRNA1, dims = pc_dims)
  scRNA1 <- FindClusters(scRNA1, resolution = cluster_resolution)
  scRNA1 <- RunUMAP(scRNA1, dims = pc_dims)
  scRNA1 <- RunTSNE(scRNA1, dims = pc_dims)
}

if (!"sample_group" %in% colnames(scRNA1@meta.data)) {
  scRNA1$sample_group <- infer_sample_group(as.character(scRNA1$orig.ident))
}

cluster_ids <- sort(unique(as.character(scRNA1$seurat_clusters)))
message("Observed Fig S1A clusters: ", length(cluster_ids), " (",
        paste(cluster_ids, collapse = ", "), ")")
if (length(cluster_ids) != 21) {
  warning("PDF Fig S1A shows 21 clusters. Current object has ",
          length(cluster_ids), ". Check resolution/source object before final figure export.")
}

p_s1a <- DimPlot(
  scRNA1,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  label.size = 4
) + ggtitle("Clusters")

ggsave("output/FigS1A_cluster_UMAP.pdf", p_s1a, width = 8, height = 7)
ggsave("output/FigS1A_cluster_UMAP.png", p_s1a, width = 8, height = 7,
       dpi = 300, bg = "white")

saveRDS(scRNA1, "output/scRNA1_preprocessed.rds")
message("Done: Fig S1A preprocessing object and cluster UMAP saved.")
