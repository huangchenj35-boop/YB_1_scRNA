## ============================================================
## Fig S3E-H — YB-1 + ABCB1/ABCC1/ABCC2/MVP 5-gene heatmap by 4 groupings
## ------------------------------------------------------------
## Source : 17_FigS3E-H_5gene_heatmap.R (5-gene group-mean z-score heatmap).
##          Layout rendered with ComplexHeatmap to match the article exactly:
##          2 square columns, gene labels on the LEFT, a thin group colour bar
##          on top, z-score legend (-2/0/2) on the right.
##          Panels: S3E sample_group, S3F group_copykat,
##                  S3G Cisplatin_pred_group, S3H Gemcitabine_pred_group.
## Output : results/FigS3{E,F,G,H}_5gene_heatmap_*.pdf / .png
## ============================================================

suppressPackageStartupMessages({
  library(Seurat); library(ComplexHeatmap); library(circlize); library(RColorBrewer); library(grid)
})
if (!exists("RESULTS")) { RESULTS <- Sys.getenv("YB1_RESULTS", unset="../results"); dir.create(RESULTS, recursive=TRUE, showWarnings=FALSE) }

genes5  <- c("YBX1","ABCB1","ABCC1","ABCC2","MVP")
col_fun <- circlize::colorRamp2(c(-2,0,2), c("#2166ac","#f7fbff","#b2182b"))
CELL    <- unit(22, "pt")   # 方格边长（版式与原图一致）

draw_one <- function(panel, obj, grp_col, grp_levels) {
  if (is.null(obj) || !inherits(obj,"Seurat")) { message("skip ",panel,": object missing"); return(invisible()) }
  if (!(grp_col %in% colnames(obj@meta.data))) { message("skip ",panel,": ",grp_col," missing"); return(invisible()) }
  DefaultAssay(obj) <- "RNA"
  gp <- genes5[genes5 %in% rownames(obj)]
  expr <- as.matrix(tryCatch(GetAssayData(obj, assay="RNA", layer="data"),
                             error=function(e) GetAssayData(obj, assay="RNA", slot="data"))[gp, , drop=FALSE])
  grp <- factor(obj@meta.data[[grp_col]], levels=grp_levels); names(grp) <- colnames(obj)
  keep <- !is.na(grp); expr <- expr[, keep, drop=FALSE]; grp <- grp[keep]
  ## 2-column group means, then row z-score (= pheatmap scale="row")
  ht <- sapply(grp_levels, function(g) rowMeans(expr[, grp==g, drop=FALSE], na.rm=TRUE))
  ht <- as.matrix(ht); colnames(ht) <- grp_levels
  z  <- t(apply(ht, 1, function(x){ s<-stats::sd(x); if(!is.finite(s)||s==0) rep(0,length(x)) else (x-mean(x))/s }))
  rownames(z) <- rownames(ht); colnames(z) <- grp_levels
  rn <- rownames(z); rn[rn=="YBX1"] <- "YB-1"; rownames(z) <- rn

  ## 与原图一致：列名斜排显示在上方（无顶部色条），行名在左，方格，z-score 图例
  ht_obj <- ComplexHeatmap::Heatmap(
    z, name = "z-score", col = col_fun,
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_names_side = "left",
    show_column_names = TRUE, column_names_side = "top", column_names_rot = 45,
    width = ncol(z)*CELL, height = nrow(z)*CELL,
    rect_gp = grid::gpar(col = "grey50", lwd = 0.8),
    heatmap_legend_param = list(at = c(-2,0,2), title = "z-score")
  )
  for (dev in c("pdf","png")) {
    fn <- file.path(RESULTS, paste0("FigS3", panel, "_5gene_heatmap_by_", grp_col, ".", dev))
    if (dev=="pdf") pdf(fn, width=3.0, height=2.8) else png(fn, width=3.0, height=2.8, units="in", res=300, bg="white")
    ComplexHeatmap::draw(ht_obj)
    dev.off()
  }
  message("done FigS3", panel)
}

## S3E sample_group needs Adjacent+Tumor -> full object scRNA1 (the epithelial
## subset is all Tumor). S3G/S3H use the aneuploid drug objects (the cells that
## carry the predicted-sensitivity grouping), matching the paper.
draw_one("E", if (exists("scRNA1")) scRNA1 else NULL,            "sample_group",          c("Tumor","Adjacent"))
draw_one("F", if (exists("scRNA_sub")) scRNA_sub else NULL,      "group_copykat",         c("aneuploid","diploid"))
draw_one("G", if (exists("scRNA_cis")) scRNA_cis else NULL,      "Cisplatin_pred_group",  c("Predicted sensitive","Predicted resistant"))
draw_one("H", if (exists("scRNA_gem")) scRNA_gem else NULL,      "Gemcitabine_pred_group",c("Predicted sensitive","Predicted resistant"))
