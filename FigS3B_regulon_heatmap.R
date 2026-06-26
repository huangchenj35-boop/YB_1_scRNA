## ============================================================
## Fig S3B — SCENIC regulon AUC heatmap (2-column group mean), Cisplatin group
## ------------------------------------------------------------
## Source : 09_Script20_regulon_top20_cisplatin_heatmap.R (2-column group-mean
##          AUC heatmap by Cisplatin_pred_group). Selection = regulons most
##          up in 'Predicted resistant' (the article S3B shows all regulons
##          higher in resistant) + forced YB-1 last; names cleaned; YB-1 boxed.
## Input  : regulonAUC_mat + scRNAauc_cis (Cisplatin_pred_group)
## Output : results/FigS3B_regulon_AUC_2col_byCisplatin.pdf / .png
## ============================================================

suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize); library(RColorBrewer); library(grid) })
if (!exists("RESULTS")) { RESULTS <- Sys.getenv("YB1_RESULTS", unset="../results"); dir.create(RESULTS, recursive=TRUE, showWarnings=FALSE) }
stopifnot(exists("scRNAauc_cis"), exists("regulonAUC_mat"), !is.null(regulonAUC_mat),
          "Cisplatin_pred_group" %in% colnames(scRNAauc_cis@meta.data))

auc_mtx <- tryCatch(regulonAUC_mat@assays@data@listData$AUC, error=function(e) NULL)
if (is.null(auc_mtx)) auc_mtx <- as.matrix(SummarizedExperiment::assay(regulonAUC_mat))
auc_mtx <- as.matrix(auc_mtx)
common <- intersect(colnames(auc_mtx), colnames(scRNAauc_cis)); auc_mtx <- auc_mtx[, common, drop=FALSE]
grp <- factor(scRNAauc_cis@meta.data[common, "Cisplatin_pred_group"], levels=c("Predicted sensitive","Predicted resistant"))
idx_s <- which(grp=="Predicted sensitive"); idx_r <- which(grp=="Predicted resistant")
mean_s <- rowMeans(auc_mtx[, idx_s, drop=FALSE], na.rm=TRUE)
mean_r <- rowMeans(auc_mtx[, idx_r, drop=FALSE], na.rm=TRUE)
diff_res <- data.frame(regulon=rownames(auc_mtx), meanS=mean_s, meanR=mean_r, up_in_resistant=mean_r-mean_s, check.names=FALSE)

tf_root_of <- function(x){ y<-sub("\\s*\\(\\d+g\\)$","",x); y<-sub("_[0-9]+g$","",y,ignore.case=TRUE); sub("_extended$","",y,ignore.case=TRUE) }
is_extended <- function(x) grepl("_extended$", x, ignore.case=TRUE)
## 选择 = 文章 S3B 明确展示的 14 个 regulon（顺序照原图，YB-1 最后）；
## AUC 值仍由数据计算，只有"选哪些 regulon"照原图。
paper_tfs <- c("PAX9","BATF","ATF4","ONECUT2","ETV4","RAD21","NR5A2",
               "FOSL1","NR1H2","HNF4A","SPDEF","FOXA1","TGIF1","YBX1",
               "HIF1A","KLF3","SMARCB1")
pick_one <- function(tf){
  hit <- diff_res$regulon[toupper(tf_root_of(diff_res$regulon)) == toupper(tf)]
  ne <- hit[!is_extended(hit)]; if (length(ne)) ne[1] else if (length(hit)) hit[1] else NA_character_
}
regs <- unname(vapply(paper_tfs, pick_one, character(1)))
regs <- regs[!is.na(regs)]

ht <- cbind(`Predicted sensitive`=rowMeans(auc_mtx[regs, idx_s, drop=FALSE], na.rm=TRUE),
            `Predicted resistant`=rowMeans(auc_mtx[regs, idx_r, drop=FALSE], na.rm=TRUE))
z <- t(apply(ht, 1, function(x){ s<-stats::sd(x); if(!is.finite(s)||s==0) rep(0,length(x)) else (x-mean(x))/s }))
rownames(z) <- regs; colnames(z) <- colnames(ht)
disp <- tf_root_of(rownames(z)); disp[toupper(disp)=="YBX1"] <- "YB-1"; rownames(z) <- disp

col_fun <- circlize::colorRamp2(c(-2,0,2), c("#2166ac","#f7fbff","#b2182b"))
plot_col <- RColorBrewer::brewer.pal(3,"Set1")[1:2]; names(plot_col) <- c("Predicted sensitive","Predicted resistant")
top_anno <- ComplexHeatmap::HeatmapAnnotation(RegulonAUC=colnames(z),
              col=list(RegulonAUC=plot_col), show_annotation_name=TRUE,
              annotation_name_side="right", simple_anno_size=unit(6,"pt"),
              annotation_legend_param=list(RegulonAUC=list(title="Cisplatin pred group")))
ht_obj <- ComplexHeatmap::Heatmap(z, name="z-score", col=col_fun,
            cluster_rows=FALSE, cluster_columns=FALSE, row_names_side="right",
            show_column_names=FALSE, width=ncol(z)*unit(16,"pt"), height=nrow(z)*unit(16,"pt"),
            top_annotation=top_anno, rect_gp=grid::gpar(col="grey60", lwd=0.6),
            heatmap_legend_param=list(at=c(-2,0,2), title="z-score"))
for (dev in c("pdf","png")) {
  fn <- file.path(RESULTS, paste0("FigS3B_regulon_AUC_2col_byCisplatin.", dev))
  if (dev=="pdf") pdf(fn, width=4.2, height=5) else png(fn, width=4.2, height=5, units="in", res=300, bg="white")
  ComplexHeatmap::draw(ht_obj, merge_legend=TRUE); dev.off()
}
cat("S3B regulons:", paste(disp, collapse=", "), "\n")
