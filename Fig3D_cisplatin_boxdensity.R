## ============================================================
## Fig 3D — YB-1 / YB-1 regulonAUC / GSVA YB-1 targets by Cisplatin group
## ------------------------------------------------------------
## Source : same raincloud (density + box) engine as Fig 1D/1G (11_Script37).
##          Grouping = Cisplatin_pred_group (Predicted sensitive/resistant).
## Input  : scRNAauc_cis (YBX1, YBX1_extended_907g, GSVA_YBX1_targets_ssgsea,
##          Cisplatin_pred_group)
## Output : results/Fig3D_<feature>_by_Cisplatin.pdf / .png
## ============================================================

suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(aplot); library(scales) })
if (!exists("RESULTS")) { RESULTS <- Sys.getenv("YB1_RESULTS", unset="../results"); dir.create(RESULTS, recursive=TRUE, showWarnings=FALSE) }
stopifnot(exists("scRNAauc_cis"), "Cisplatin_pred_group" %in% colnames(scRNAauc_cis@meta.data))

## GSVA_YBX1_targets_ssgsea lives in scRNAauc (same 12290 cells); copy by barcode if absent
if (!("GSVA_YBX1_targets_ssgsea" %in% colnames(scRNAauc_cis@meta.data)) &&
    exists("scRNAauc") && "GSVA_YBX1_targets_ssgsea" %in% colnames(scRNAauc@meta.data)) {
  g <- scRNAauc@meta.data[["GSVA_YBX1_targets_ssgsea"]]; names(g) <- colnames(scRNAauc)
  scRNAauc_cis$GSVA_YBX1_targets_ssgsea <- g[colnames(scRNAauc_cis)]
}

lev <- c("Predicted sensitive","Predicted resistant")
pal <- c("Predicted sensitive"="#5AAFD6", "Predicted resistant"="#F29C64")  # 敏感蓝/耐药橙

raincloud <- function(obj, feature, title, fout){
  v <- if (feature %in% colnames(obj@meta.data)) obj@meta.data[[feature]] else FetchData(obj, vars=feature)[,1]
  grp <- factor(as.character(obj$Cisplatin_pred_group), levels=lev)
  df <- data.frame(Subtype=grp, AUC=as.numeric(v)); df <- df[is.finite(df$AUC) & !is.na(df$Subtype), ]
  p.val <- tryCatch(wilcox.test(AUC ~ Subtype, data=df)$p.value, error=function(e) NA_real_)
  p.lab <- if (is.finite(p.val)) paste0("P", ifelse(p.val<0.001," < 0.001", paste0(" = ", round(p.val,3)))) else "P = NA"
  fill_vals <- pal[lev]; alpha_vals <- alpha(fill_vals,0.7); names(alpha_vals) <- lev
  p_top <- ggplot(df, aes(x=AUC, color=Subtype, fill=Subtype)) + geom_density() +
    scale_color_manual(values=alpha_vals, na.translate=FALSE) + scale_fill_manual(values=alpha_vals, na.translate=FALSE) +
    labs(title=title, x=NULL, y=NULL) +
    annotate("text", x=Inf, y=Inf, label="Predicted resistant", color=pal[["Predicted resistant"]], hjust=1.05, vjust=1.6, size=3.8) +
    annotate("text", x=Inf, y=Inf, label="Predicted sensitive", color=pal[["Predicted sensitive"]], hjust=1.05, vjust=3.0, size=3.8) +
    theme_classic() + theme(legend.position="none", plot.title=element_text(hjust=0.5,size=14),
      axis.text.x=element_text(size=12,color="black"), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
      axis.line.y=element_blank(), panel.grid=element_blank(), plot.margin=margin(4,40,2,4)) + geom_rug()
  p_bot <- ggplot(df, aes(Subtype, AUC, fill=Subtype)) + geom_boxplot(aes(col=Subtype)) +
    scale_fill_manual(values=fill_vals, na.translate=FALSE) + scale_color_manual(values=fill_vals, na.translate=FALSE) +
    xlab(NULL) + ylab(NULL) + theme_void() +
    theme(legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank(),
          panel.grid=element_blank(), plot.margin=margin(2,40,2,4)) +
    annotate("text", x=1.5, hjust=0.5, y=max(df$AUC,na.rm=TRUE)*0.97, size=3.3, angle=270, fontface="bold", label=p.lab) +
    coord_flip(clip="off")
  dat <- ggplot_build(p_bot)$data[[1]]
  if (!is.null(dat) && nrow(dat)) p_bot <- p_bot + geom_segment(data=dat, aes(x=xmin,xend=xmax,y=middle,yend=middle), color="white", inherit.aes=FALSE)
  p_all <- p_top %>% insert_bottom(p_bot, height=0.45)
  ggsave(file.path(RESULTS, paste0(fout,".pdf")), p_all, width=3.4, height=3.4)
  ggsave(file.path(RESULTS, paste0(fout,".png")), p_all, width=3.4, height=3.4, dpi=300, bg="white")
  message("done ", fout)
}

raincloud(scRNAauc_cis, "YBX1",                     "YB-1",            "Fig3D_YBX1_by_Cisplatin")
raincloud(scRNAauc_cis, "YBX1_extended_907g",       "YB-1 regulonAUC", "Fig3D_YBX1regulonAUC_by_Cisplatin")
raincloud(scRNAauc_cis, "GSVA_YBX1_targets_ssgsea", "GSVA YB-1 targets","Fig3D_GSVA_by_Cisplatin")
