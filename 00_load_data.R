## ============================================================
## 00_load_data.R  —  Code Ocean data loader
## ------------------------------------------------------------
## Per the agreed rule (Choice A): this file ONLY sets paths and
## loads the already-processed objects into the variable names the
## original scripts assume in memory. It does NOT change any
## analysis or plotting logic of the figure scripts.
##
## Usage (inside R, before sourcing a figure script):
##     source("00_load_data.R")
##     source("Fig1B_sample_origin_UMAP.R")
##
## Object -> figure map (see RUN_ORDER.md):
##   scRNA1    (full 33990)      -> Fig 1A, 1B, 1C, 1D, S1C
##   scRNA_sub (epithelial 12290)-> Fig 1E, 1F, 1G
##   scRNA_cis / scRNA_gem       -> Fig 3D / 3E (pending call)
##
## On Code Ocean the data are mounted at /data; locally they sit in
## ../data. Override with environment variable YB1_DATA if needed.
## ============================================================

## Attach the libraries the rougher scripts assume are already loaded
## (some original scripts call DimPlot/ggsave without their own library() line).
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

DATA <- Sys.getenv("YB1_DATA", unset = "../data")
if (!dir.exists(DATA)) stop("Data directory not found: ", DATA,
                            "  (set env var YB1_DATA to the data folder)")

## Output directory: all figures are written here, named in article order
## (Fig1B_..., Fig1G_..., FigS1C_...). Code Ocean mounts results at /results.
RESULTS <- Sys.getenv("YB1_RESULTS", unset = "../results")
dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)

## ---- 1. Main curated single-cell object  ->  scRNA1 -----------------
## meta.data contains: curate_v1 (13 cell types), sample_group
## (Adjacent/Tumor), group_copykat (aneuploid/diploid), seurat_clusters,
## curate_malig.  Assay "RNA" contains gene YBX1.  Reductions: pca, umap.
main_rdata <- file.path(DATA,
  "integrated_FindClusters_umap_singleR_0.5_curate_cnv_malignant_2.Rdata")
if (file.exists(main_rdata)) {
  load(main_rdata)                       # restores object named 'scRNA1'
  if (!"sample_type" %in% colnames(scRNA1@meta.data))
    scRNA1$sample_type <- scRNA1$sample_group   # some scripts read sample_type
  message("Loaded scRNA1: ", ncol(scRNA1), " cells, ",
          nrow(scRNA1), " genes")
} else {
  warning("Main object missing: ", main_rdata)
}

## ---- 2. Aneuploid subset + drug-sensitivity groups ------------------
## scRNA_cis: meta has Cisplatin_pred_group, IC50uM_Cisplatin, lnIC50_Cisplatin
## scRNA_gem: meta has Gemcitabine_pred_group, IC50uM_Gemcitabine, lnIC50_Gemcitabine
cis_rds <- file.path(DATA, "scRNA_aneuploid_with_Cisplatin_groups.rds")
gem_rds <- file.path(DATA, "scRNA_aneuploid_with_Gemcitabine_groups.rds")
scRNA_cis <- if (file.exists(cis_rds)) readRDS(cis_rds) else NULL
scRNA_gem <- if (file.exists(gem_rds)) readRDS(gem_rds) else NULL

## ---- 2b. Epithelial (cholangiocyte + hepatocyte) subset, own UMAP ------
## scRNA_sub: 12290 epithelial cells (Cholangiocytes 12250 + Hepatocytes 40),
## meta col group_copykat (aneuploid 7040 / diploid 5250), gene YBX1,
## reductions pca/umap/tsne = the subset's own re-embedding.
## This is the embedding used in the paper for Fig 1C, Fig 1E, Fig 1F.
sub_rdata <- file.path(DATA, "cholangiocyte_subset_FindClusters_sct_umap.Rdata")
scRNA_sub <- NULL
if (file.exists(sub_rdata)) {
  se <- new.env(); load(sub_rdata, envir = se)
  scRNA_sub <- get(ls(se)[1], envir = se)   # object inside is named 'scRNA1'
  if ("RNA" %in% names(scRNA_sub@assays)) DefaultAssay(scRNA_sub) <- "RNA"
}

## ---- 3. SCENIC / GSVA objects (epithelial subset, 12290 cells) ------
## Unlock Fig 3A/3B/3C, Fig S3A/S3B/S3C/S3D and Fig S3E-H.
##  scRNAauc      : meta YBX1_extended_907g (regulon AUC) + GSVA_YBX1_targets_ssgsea
##                  + META assay "YBX1-extended-907g"; reductions umap/tsne
##                  (scRNAauc 自带 GSVA_YBX1_targets_ssgsea，供 Fig 3C/3D/3E)
##  scRNAauc_cis  : adds lnIC50_Cisplatin_1005 + Cisplatin_pred_group to scRNAauc
##                  (Fig S3C density, Fig 3D, Fig S3A/S3B by cisplatin group)
auc_rds   <- file.path(DATA, "scRNAauc.rds")
auccis_rds<- file.path(DATA, "scRNAauc_with_CisplatinGroup_Cisplatin.rds")
scRNAauc    <- if (file.exists(auc_rds))   readRDS(auc_rds)    else NULL
scRNAauc_cis<- if (file.exists(auccis_rds))readRDS(auccis_rds) else NULL

## SCENIC regulon AUC / binary matrices (for Fig 3A heatmap & Fig S3A volcano)
regAUC_rds <- file.path(DATA, "3.4_regulonAUC.Rds")
regBIN_rds <- file.path(DATA, "4.1_binaryRegulonActivity.Rds")
regulonAUC_mat <- if (file.exists(regAUC_rds)) readRDS(regAUC_rds) else NULL
binaryReg_mat  <- if (file.exists(regBIN_rds)) readRDS(regBIN_rds) else NULL

message("Loader done. Available objects: ",
        paste(c("scRNA1",
                if(!is.null(scRNA_cis))    "scRNA_cis",
                if(!is.null(scRNA_gem))    "scRNA_gem",
                if(!is.null(scRNA_sub))    "scRNA_sub",
                if(!is.null(scRNAauc))     "scRNAauc",
                if(!is.null(scRNAauc_cis)) "scRNAauc_cis",
                if(!is.null(regulonAUC_mat)) "regulonAUC_mat",
                if(!is.null(binaryReg_mat))  "binaryReg_mat"), collapse = ", "))
