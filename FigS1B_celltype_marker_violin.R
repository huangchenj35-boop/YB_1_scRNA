## ============================================================
## Fig S1B — Cell-type marker stacked violin (full object)
## ------------------------------------------------------------
## Source : author's marker-violin template (provided verbatim).
##          Adapted only: object name scRNA1_iCCA -> scRNA1; Seurat-5 data
##          accessor; full article gene list + column order; output path.
##          Plotting (p1) is the author's template, unchanged.
## Input  : scRNA1 (loaded by 00_load_data.R; meta col curate_v1, RNA assay)
## Output : results/FigS1B_celltype_marker_violin.pdf / .png
##          Column colours = ggplot default hue in Fig1A order (matches Fig1A).
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)        # 原脚本 library(tidyverse) 中实际用到的部分
  library(ggplot2)
  library(ggbeeswarm)
  library(scales)
  library(reshape2)
})

if (!exists("RESULTS")) {
  RESULTS <- Sys.getenv("YB1_RESULTS", unset = "../results")
  dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
}

stopifnot(exists("scRNA1"), "curate_v1" %in% colnames(scRNA1@meta.data))

## ---- article Fig S1B marker gene list (top -> bottom) ----
features = c("CD19","MS4A1","CD79A","CD79B","BANK1","LAG3","PDCD1","SOX9","KRT19",
             "EPCAM","CLEC9A","XCR1","BATF3","IRF8","THBD","ECSCR","VWF","ESAM",
             "COL1A2","COL6A1","MYL9","COL1A1","MMP2","MCAM","IL6","PDGFA","MYLK",
             "COL3A1","ALB","CYP3A4","APOA1","HNF4A","FCGR2A","CD86","CD14","FCGR3A",
             "S100A8","S100A9","CD68","CSF1R","IFNG","KLRC1","KLRD1","KLRB1","IGJ",
             "SLAMF7","LILRA4","IL3RA","TCF4","CLEC4C","IRF7","CD3D","CD3E")
features = features[features %in% rownames(scRNA1[["RNA"]])]

expr_data = GetAssayData(scRNA1, assay = "RNA", layer = "data")[features, ]   # Seurat5: was scRNA1_iCCA@assays$RNA@data[features,]
meta_data = scRNA1@meta.data
meta_data$Barcode = rownames(meta_data)
expr_data_df = as.data.frame(t(expr_data))
expr_data_df$Barcode = row.names(expr_data_df)
expr_data_df = merge(expr_data_df, meta_data, by = "Barcode")
selected_columns = c(1:(length(features) + 1), which(colnames(expr_data_df) == "curate_v1"))
expr_data_df = expr_data_df[ ,selected_columns]
expr_data_long = reshape2::melt(expr_data_df, id.vars = c("Barcode", "curate_v1"),
                                variable.names = "Gene", value.name = "Expression")
colnames(expr_data_long)[colnames(expr_data_long) == "Barcode"] = "CB"
colnames(expr_data_long)[colnames(expr_data_long) == "curate_v1"] = "celltype"
colnames(expr_data_long)[colnames(expr_data_long) == "variable"] = "gene"
colnames(expr_data_long)[colnames(expr_data_long) == "Expression"] = "exp"
vln.df = expr_data_long

## display label + column order to match the article (CD8+ Tex cell -> CD8+ T cell)
vln.df$celltype = as.character(vln.df$celltype)
vln.df$celltype[vln.df$celltype == "CD8+ Tex cell"] = "CD8+ T cell"
celltype_order = c("B cell","CD8+ T cell","Cholangiocytes","Dendritic cells",
                   "Endothelial cell","Fibroblast","Hepatocytes","moDCs","Mono/Macro",
                   "NK cells","Plasma cell","Plasmacytoid dendritic cells","T cells")
vln.df$celltype = factor(vln.df$celltype, levels = celltype_order)
vln.df$gene = factor(vln.df$gene, levels = features)

p1 <- vln.df%>%ggplot(aes(celltype,exp),color=factor(celltype))+
  geom_violin(aes(fill=celltype),scale = "width")+
  facet_grid(gene~.,scales = "free_y")+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),

    axis.title.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = NULL,color = "black",size = 14),
    axis.title.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.text.y.left = element_blank(),

    legend.position = "none",

    panel.spacing.y = unit(0, "cm"),
    strip.text.y = element_text(angle=0,size = 14,hjust = 0),
    strip.background.y = element_blank(),
    plot.margin = unit(c(0.2, 0.2, 1.6, 0.2), "cm")   # 底部留白，避免细胞类型长标签被截断
  )
p1
ggsave(file.path(RESULTS, "FigS1B_celltype_marker_violin.pdf"), p1, width = 12, height = 34, units = "cm", limitsize = FALSE)
ggsave(file.path(RESULTS, "FigS1B_celltype_marker_violin.png"), p1, width = 12, height = 34, units = "cm", dpi = 300, bg = "white", limitsize = FALSE)
