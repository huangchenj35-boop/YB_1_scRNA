## ============================================================
## Fig S1C — Cell-type composition (Adjacent vs Tumor + per sample)
## ------------------------------------------------------------
## Input : scRNA1   (loaded by 00_load_data.R; meta cols curate_v1, orig.ident)
## Output: results/FigS1C_Stacked_Proportion_AT.pdf / .png            (Fig S1C)
##         results/FigS1C_Stacked_Proportion_by_Sample.pdf / .png     (per sample)
##         results/FigS1C_CellType_Proportion_by_Sample.csv           (table)
## Note  : plotting parameters unchanged; only header/guard added.
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(colorspace)   # 用于把既定调色板整体“变浅”，不改颜色体系
})

# allow standalone use if 00_load_data.R was not sourced
if (!exists("RESULTS")) {
  RESULTS <- Sys.getenv("YB1_RESULTS", unset = "../results")
  dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
}

## ========== 0) 准备对象 & 列名检查 ==========
stopifnot(exists("scRNA1"), inherits(scRNA1, "Seurat"))
need_cols <- c("curate_v1","orig.ident")
miss <- setdiff(need_cols, colnames(scRNA1@meta.data))
if (length(miss)) stop("meta.data 缺少列：", paste(miss, collapse = ", "))

## 细胞顺序（严格按此显示；显示名与文章一致："CD8+ T cell"）
cell_order <- c(
  "B cell","CD8+ T cell","Cholangiocytes","Dendritic cells","Endothelial cell",
  "Fibroblast","Hepatocytes","moDCs","Mono/Macro","NK cells",
  "Plasma cell","Plasmacytoid dendritic cells","T cells"
)

## 你给的配色
## (restored from the author's commented line; script references zzm8colors below)
zzm8colors  <- c('#d4de9c','#94c58f','#86c7b4','#9cd2ed','#a992c0',
                 '#ea9994','#f2c396','#bb82b1')
zzm60colors <- c('#4b6aa8','#3ca0cf','#c376a7','#ad98c3','#cea5c7',
                 '#53738c','#a5a9b0','#a78982','#92699e',
                 '#d69971','#df5734','#6c408e','#ac6894','#d4c2db',
                 '#537eb7','#83ab8e','#ece399','#405993','#cc7f73',
                 '#b95055','#d5bb72','#bc9a7f','#e0cfda','#d8a0c0',
                 '#e6b884','#b05545','#d69a55','#64a776','#cbdaa9',
                 '#efd2c9','#da6f6d','#ebb1a4','#a44e89','#a9c2cb',
                 '#b85292','#6d6fa0','#8d689d','#c8c7e1','#d25774',
                 '#c49abc','#927c9a','#3674a2','#9f8d89','#72567a',
                 '#63a3b8','#c4daec','#61bada','#b7deea','#e29eaf',
                 '#4490c4','#e6e2a3','#de8b36','#c4612f','#9a70a8',
                 '#76a2be','#408444','#c6adb0','#9d3b62','#2d3462')

## ========== 1) 锁定细胞类型顺序 & 只保留这 13 类 ==========
# 显示名对齐文章：对象里的 "CD8+ Tex cell" 在图中显示为 "CD8+ T cell"
cv_disp <- as.character(scRNA1$curate_v1)
cv_disp[cv_disp == "CD8+ Tex cell"] <- "CD8+ T cell"
scRNA1$curate_v1 <- factor(cv_disp, levels = cell_order)
cells_keep <- rownames(scRNA1@meta.data)[!is.na(scRNA1$curate_v1)]
scRNA1_use <- subset(scRNA1, cells = cells_keep)

## 固定样本顺序（给定那 8 个样本，缺失的自动忽略）
sample_order_given <- c(
  "GSM4116579_ICC_18_Adjacent",
  "GSM4116580_ICC_18_Tumor",
  "GSM4116581_ICC_20_Tumor",
  "GSM4116582_ICC_23_Adjacent",
  "GSM4116583_ICC_23_Tumor",
  "GSM4116584_ICC_24_Tumor1",
  "GSM4116585_ICC_24_Tumor2",
  "GSM4116586_ICC_25_Adjacent"
)
present_samples <- intersect(sample_order_given, unique(scRNA1_use$orig.ident))
if (length(present_samples) == 0) present_samples <- unique(scRNA1_use$orig.ident)
scRNA1_use$orig.ident <- factor(as.character(scRNA1_use$orig.ident), levels = present_samples)

## ========== 2) 计数 & 百分比 ==========
meta_df <- scRNA1_use@meta.data %>%
  dplyr::select(orig.ident, curate_v1)

count_df <- meta_df %>%
  count(orig.ident, curate_v1, name = "n_cells") %>%
  complete(orig.ident, curate_v1, fill = list(n_cells = 0)) %>%
  group_by(orig.ident) %>%
  mutate(total_cells = sum(n_cells),
         percent = 100 * n_cells / total_cells) %>%
  ungroup()

## 输出表格
write.csv(count_df, file.path(RESULTS, "FigS1C_CellType_Proportion_by_Sample.csv"), row.names = FALSE)

## ========== 3) 配色：直接采用文章原图配色（从 PDF 图例逐格取色） ==========
pal_use <- c(
  "B cell"                       = "#A2C8DC",
  "CD8+ T cell"                  = "#2A79AC",
  "Cholangiocytes"               = "#93C58D",
  "Dendritic cells"              = "#50A643",
  "Endothelial cell"             = "#B49973",
  "Fibroblast"                   = "#DF504E",
  "Hepatocytes"                  = "#E26A45",
  "moDCs"                        = "#F19E42",
  "Mono/Macro"                   = "#E18947",
  "NK cells"                     = "#AC90C0",
  "Plasma cell"                  = "#7B5A92",
  "Plasmacytoid dendritic cells" = "#F0EB8E",
  "T cells"                      = "#AA5629"
)

## ========== 4) 每个样本的细胞组成（堆叠比例柱状图） ==========
p_stack <- ggplot(count_df, aes(x = orig.ident, y = percent, fill = curate_v1)) +
  geom_col(width = 0.9, color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = pal_use, name = "Cell type") +
  labs(x = NULL, y = "Proportion (%)", title = "Cell-type composition per sample") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10)
  )

ggsave(file.path(RESULTS, "FigS1C_Stacked_Proportion_by_Sample.pdf"), p_stack, width = 10, height = 4.8)   # 英寸
ggsave(file.path(RESULTS, "FigS1C_Stacked_Proportion_by_Sample.png"), p_stack, width = 10, height = 4.8, dpi = 300)

## ========== 5) （可选）Adjacent vs Tumor 合并 ==========
AT_from_sample <- function(x){
  ifelse(grepl("Adjacent|Normal|Para", x, ignore.case = TRUE), "Adjacent", "Tumor")
}
count_AT <- count_df %>%
  mutate(AT = factor(AT_from_sample(as.character(orig.ident)), levels = c("Adjacent","Tumor"))) %>%
  group_by(AT, curate_v1) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop_last") %>%
  group_by(AT) %>%
  mutate(total_cells = sum(n_cells),
         ratio = n_cells / total_cells) %>%
  ungroup()

# 与文章原图一致：Y 轴为 Ratio(0–1)、图例名 "Cluster"、无标题
p_stack_AT <- ggplot(count_AT, aes(x = AT, y = ratio, fill = curate_v1)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = pal_use, name = "Cluster") +
  labs(x = NULL, y = "Ratio") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10)
  )

ggsave(file.path(RESULTS, "FigS1C_Stacked_Proportion_AT.pdf"), p_stack_AT, width = 5.6, height = 4.4)
ggsave(file.path(RESULTS, "FigS1C_Stacked_Proportion_AT.png"), p_stack_AT, width = 5.6, height = 4.4, dpi = 300)

cat("✅ 已完成：\n",
    "1) 样本层面的 13 类细胞百分比统计（CellType_Proportion_by_Sample.csv）。\n",
    "2) 堆叠比例图：Stacked_Proportion_by_Sample.pdf/.png（统一浅色化 zzm 调色）。\n",
    "3) Adjacent vs Tumor 合并对比：Stacked_Proportion_AT.pdf/.png（同一套配色）。\n", sep = "")
