## ============================================================
## Fig 3A — SCENIC regulon AUC per-cell heatmap, aneuploid vs diploid
## ------------------------------------------------------------
## Source : 06_Script18_SCENIC_AUC_BIN_heatmap.R (the AUC heatmap block,
##          show_colnames=FALSE -> per-cell columns; ported verbatim).
##          Adapted: AUCmatrix from data/3.4_regulonAUC.Rds (regulon x cell,
##          names cleaned), scRNA <- scRNAauc, output -> results/.
## Input  : regulonAUC_mat + scRNAauc (group_copykat)
## Output : results/Fig3A_AUC_heatmap_selected_regulons_YBX1_Top20.pdf / .png
## ============================================================

suppressPackageStartupMessages({ library(Seurat); library(pheatmap); library(RColorBrewer) })
if (!exists("RESULTS")) { RESULTS <- Sys.getenv("YB1_RESULTS", unset="../results"); dir.create(RESULTS, recursive=TRUE, showWarnings=FALSE) }
stopifnot(exists("scRNAauc"), exists("regulonAUC_mat"), !is.null(regulonAUC_mat))

## AUCmatrix = regulon x cell, names cleaned exactly as 06_Script18 (" (xxxg)" -> "_xxxg")
AUCmatrix <- tryCatch(regulonAUC_mat@assays@data@listData$AUC, error=function(e) NULL)
if (is.null(AUCmatrix)) AUCmatrix <- as.matrix(SummarizedExperiment::assay(regulonAUC_mat))
AUCmatrix <- as.matrix(AUCmatrix)
rn <- rownames(AUCmatrix); rn <- gsub(' \\(','_',rn); rn <- gsub('\\)','',rn); rownames(AUCmatrix) <- rn
scRNA <- scRNAauc
common <- intersect(colnames(AUCmatrix), colnames(scRNA)); AUCmatrix <- AUCmatrix[, common, drop=FALSE]

## ===== 06_Script18 selection (unchanged) =====
tf_root_of <- function(x){ y<-sub("\\s*\\(\\d+g\\)$","",x); y<-sub("_[0-9]+g$","",y,ignore.case=TRUE); sub("_extended$","",y,ignore.case=TRUE) }
is_extended <- function(x) grepl("_extended$", x, ignore.case=TRUE)
dedup_keep_one <- function(reg_names){ if(!length(reg_names)) return(character(0))
  df<-data.frame(reg=reg_names, TF_root=tf_root_of(reg_names), ext=is_extended(reg_names), stringsAsFactors=FALSE)
  df<-df[order(df$TF_root, df$ext), ]; df$reg[!duplicated(df$TF_root)] }
seed_regs <- rownames(AUCmatrix)[80:100]
regs_available <- rownames(AUCmatrix)
regs_dedup <- dedup_keep_one(intersect(seed_regs, regs_available))
ybx1_all <- grep("^YBX1(\\b|_|\\s|\\()", rownames(AUCmatrix), value=TRUE)
ybx1_pick <- if(length(ybx1_all)){ if(any(!is_extended(ybx1_all))) ybx1_all[!is_extended(ybx1_all)][1] else ybx1_all[1] } else NA_character_
all_dedup <- dedup_keep_one(regs_available)
prior <- unique(na.omit(c(ybx1_pick, regs_dedup))); fill <- setdiff(all_dedup, prior)
regs20 <- unique(c(prior, fill))[1:min(20, length(unique(c(prior, fill))))]
regs20 <- dedup_keep_one(regs20)

## ===== per-cell AUC heatmap (unchanged) =====
myAUCmatrix <- AUCmatrix[regs20, , drop=FALSE]
grp <- scRNA@meta.data[colnames(myAUCmatrix), "group_copykat"]
grp <- factor(grp, levels = intersect(c("aneuploid","diploid"), unique(as.character(grp))))
keep <- !is.na(grp); myAUCmatrix <- myAUCmatrix[, keep, drop=FALSE]; grp <- droplevels(grp[keep])
## order cells aneuploid then diploid (so the two blocks are contiguous)
ord_cells <- order(grp)
myAUCmatrix <- myAUCmatrix[, ord_cells, drop=FALSE]; grp <- grp[ord_cells]
anno_col <- data.frame(group_copykat = grp); rownames(anno_col) <- colnames(myAUCmatrix)
plot_col <- brewer.pal(3,"Set1")[1:2]; names(plot_col) <- c("aneuploid","diploid")
annotation_colors <- list(group_copykat = plot_col)
row_z <- function(m){ out<-t(apply(m,1,function(x){mu<-mean(x,na.rm=TRUE);sdv<-stats::sd(x,na.rm=TRUE); if(!is.finite(sdv)||sdv==0) rep(0,length(x)) else (x-mu)/sdv})); rownames(out)<-rownames(m); colnames(out)<-colnames(m); out }
myAUC_z <- row_z(myAUCmatrix)
## clean row labels to match the paper (strip _extended/_Ng; YBX1 -> YB-1)
disp <- tf_root_of(rownames(myAUC_z)); disp[toupper(disp) == "YBX1"] <- "YB-1"
rownames(myAUC_z) <- disp
bk <- c(seq(-3,-0.1,by=0.01), seq(0,3,by=0.01))
heat_colors <- c(colorRampPalette(c("#2166ac","#f7fbff"))(length(bk)/2), colorRampPalette(c("#f7fbff","#b2182b"))(length(bk)/2))
for (dev in c("pdf","png")) {
  fn <- file.path(RESULTS, paste0("Fig3A_AUC_heatmap_selected_regulons_YBX1_Top20.", dev))
  if (dev=="pdf") pdf(fn, width=8, height=6) else png(fn, width=8, height=6, units="in", res=300, bg="white")
  pheatmap(myAUC_z, show_colnames=FALSE, cluster_cols=FALSE, cluster_rows=FALSE,
           annotation_col=anno_col, annotation_colors=annotation_colors, scale="none",
           color=heat_colors, breaks=bk, border_color=NA, legend_breaks=seq(-2,2,2), use_raster=TRUE)
  dev.off()
}
cat("regs20:", paste(regs20, collapse=", "), "\n")
