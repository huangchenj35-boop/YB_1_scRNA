## ============================================================
## Fig S3A — SCENIC regulon volcano, Cisplatin sensitive vs resistant
## ------------------------------------------------------------
## Source : 14_Script15_FigS3A_volcano.R (volcano style ported verbatim).
##          Adapted only: grouping column group_copykat -> Cisplatin_pred_group
##          (the paper S3A contrast), auc_mtx from data/3.4_regulonAUC.Rds,
##          scRNA <- scRNAauc_cis, output -> results/.
## Input  : regulonAUC_mat + scRNAauc_cis (meta Cisplatin_pred_group)
## Output : results/FigS3A_SCENIC_Regulon_Volcano_Cisplatin.pdf / .png
## ============================================================

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(ggrepel)
})
if (!exists("RESULTS")) { RESULTS <- Sys.getenv("YB1_RESULTS", unset="../results"); dir.create(RESULTS, recursive=TRUE, showWarnings=FALSE) }
stopifnot(exists("scRNAauc_cis"), exists("regulonAUC_mat"), !is.null(regulonAUC_mat),
          "Cisplatin_pred_group" %in% colnames(scRNAauc_cis@meta.data))

auc_mtx <- tryCatch(regulonAUC_mat@assays@data@listData$AUC, error=function(e) NULL)
if (is.null(auc_mtx)) auc_mtx <- as.matrix(SummarizedExperiment::assay(regulonAUC_mat))
auc_mtx <- as.matrix(auc_mtx)
scRNA <- scRNAauc_cis; suffix <- ""
common <- intersect(colnames(auc_mtx), colnames(scRNA)); auc_mtx <- auc_mtx[, common, drop=FALSE]

## ---- group: Cisplatin sensitive vs resistant ----
group <- scRNA@meta.data[colnames(auc_mtx), "Cisplatin_pred_group"]
group <- factor(group, levels = c("Predicted resistant", "Predicted sensitive"))
keep <- !is.na(group); auc_mtx <- auc_mtx[, keep, drop=FALSE]; group <- droplevels(group[keep])
idx_r <- which(group == "Predicted resistant")    # a
idx_s <- which(group == "Predicted sensitive")     # d
mean_r <- rowMeans(auc_mtx[, idx_r, drop=FALSE], na.rm=TRUE)
mean_s <- rowMeans(auc_mtx[, idx_s, drop=FALSE], na.rm=TRUE)
meanDiff <- mean_r - mean_s   # >0 = up in Resistant (right); <0 = up in Sensitive (left)
pvals <- apply(auc_mtx, 1, function(v){ x<-v[idx_r]; y<-v[idx_s]
  if (all(is.na(x))||all(is.na(y))) return(NA_real_)
  tryCatch(stats::wilcox.test(x,y,exact=FALSE)$p.value, error=function(e) NA_real_) })
FDR <- p.adjust(pvals, method="BH")
df_pval <- data.frame(Regulon=rownames(auc_mtx), meanDiff=meanDiff[rownames(auc_mtx)], FDR=FDR, check.names=FALSE)
write.csv(df_pval, file.path(RESULTS, "FigS3A_SCENIC_RegulonAUC_Diff_Cisplatin.csv"), row.names=FALSE)

## ===== 14_Script15 volcano (unchanged style) =====
x <- df_pval[, c("Regulon","meanDiff","FDR")]; colnames(x) <- c("X","logFC","P.Value"); x$label <- x$X; rownames(x) <- x$X
x$TF_root <- sub("\\s*\\(\\d+g\\)$","",x$X); x$is_extended <- grepl("_extended$",x$X,ignore.case=TRUE)
x$TF_root <- sub("_extended$","",x$TF_root,ignore.case=TRUE)
x_ord <- x[order(x$P.Value, -abs(x$logFC), x$is_extended), ]; x_dedup <- x_ord[!duplicated(x_ord$TF_root), ]
ord_dedup <- order(x_dedup$P.Value, -abs(x_dedup$logFC))
top10_names <- head(rownames(x_dedup)[ord_dedup], 10)
ybx1_row <- rownames(subset(x_dedup, TF_root=="YBX1"))[1]
sel_names <- unique(na.omit(c(top10_names, ybx1_row)))
selectgenes <- x_dedup[sel_names, , drop=FALSE]; selectgenes$gsym <- rownames(selectgenes)
selectgenes$pathway <- ifelse(sub("_extended$","",sub("\\s*\\(\\d+g\\)$","",rownames(selectgenes)),ignore.case=TRUE)=="YBX1","YBX1",
                              ifelse(selectgenes$logFC>0,"Up","Down"))
logFCcut<-0.001; pvalCut<-0.05; logFCcut2<-2.5; logFCcut3<-5; pvalCut2<-1e-4; pvalCut3<-1e-5
n1<-nrow(x_dedup); cols<-rep("grey",n1); names(cols)<-rownames(x_dedup)
cols[x_dedup$P.Value<pvalCut & x_dedup$logFC> logFCcut ]<-"#FB9A99"
cols[x_dedup$P.Value<pvalCut2& x_dedup$logFC> logFCcut2]<-"#ED4F4F"
cols[x_dedup$P.Value<pvalCut & x_dedup$logFC< -logFCcut ]<-"#B2DF8A"
cols[x_dedup$P.Value<pvalCut2& x_dedup$logFC< -logFCcut2]<-"#329E3F"
x_dedup$color_transparent<-adjustcolor(cols,alpha.f=0.5)
size<-rep(1,n1); size[x_dedup$P.Value<pvalCut & abs(x_dedup$logFC)>logFCcut ]<-2
size[x_dedup$P.Value<pvalCut2& abs(x_dedup$logFC)>logFCcut2]<-4; size[x_dedup$P.Value<pvalCut3& abs(x_dedup$logFC)>logFCcut3]<-6
ymax<-max(-log10(pmax(x_dedup$P.Value,.Machine$double.xmin)))*1.1
mycol<-c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
lab_cols<-setNames(mycol[seq_along(unique(selectgenes$pathway))],unique(selectgenes$pathway))
p1<-ggplot(x_dedup, aes(logFC,-log10(P.Value),label=label))+
  geom_point(alpha=0.6,size=size,colour=x_dedup$color_transparent)+
  labs(x=bquote(~Log[2]~"(fold change of AUC: Sensitive" %<-% phantom() %->% "Resistant)"),
       y=bquote(~-Log[10]~italic("FDR")), title="")+
  ylim(c(0,ymax))+
  scale_x_continuous(breaks=c(-0.05,-0.001,0,0.001,0.05),labels=c(-0.05,-0.001,0,0.001,0.05),limits=c(-0.05,0.05))+
  geom_vline(xintercept=c(-logFCcut,logFCcut),color="grey40",linetype="longdash",linewidth=0.5)+
  geom_hline(yintercept=-log10(pvalCut),color="grey40",linetype="longdash",linewidth=0.5)+
  theme_bw(base_size=12)+theme(panel.grid=element_blank())
p2<-p1+geom_point(data=selectgenes,alpha=1,size=4.6,shape=1,stroke=1,color="black")+
  scale_color_manual(values=lab_cols)+
  geom_text_repel(data=selectgenes,aes(color=pathway),show.legend=FALSE,size=5,
                  box.padding=unit(0.35,"lines"),point.padding=unit(0.3,"lines"),max.overlaps=1000)+
  guides(color=guide_legend(title=NULL))
ggsave(file.path(RESULTS,"FigS3A_SCENIC_Regulon_Volcano_Cisplatin.pdf"), p2, width=7, height=6)
ggsave(file.path(RESULTS,"FigS3A_SCENIC_Regulon_Volcano_Cisplatin.png"), p2, width=7, height=6, dpi=300, bg="white")
