## ============================================================
## Sample-Type Visualization
## Figures: Fig 1B (sample-type UMAP), Fig 1D (YBX1 rain-cloud by sample type),
##          Fig S1C (stacked cell-type proportion bars)
## Dataset: GSE138709 (5 tumor + 3 adjacent iCCA samples)
## ============================================================
## Prerequisites:
##   - scRNA1: Seurat object with curate_v1 annotation and orig.ident
##             produced by 01_Fig1A_FigS1B_cell_annotation.R
## Output:
##   - sample_type_UMAP.pdf/.png            (Fig 1B)
##   - BoxDensity_YBX1expr_by_sample_type_withDots.pdf/.png (Fig 1D)
##   - Stacked_Proportion_by_Sample.pdf/.png (per-sample composition)
##   - Stacked_Proportion_AT.pdf/.png        (Fig S1C: Adjacent vs Tumor)
##   - CellType_Proportion_by_Sample.csv     (supplementary table)
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(colorspace)
})

stopifnot(exists("scRNA1"), inherits(scRNA1, "Seurat"))

# ------------------------------------------------------------------
# 1. Fig 1B – UMAP coloured by sample type (Tumor / Adjacent)
# ------------------------------------------------------------------
pal_contrast <- c("Adjacent" = "#0072B2",   # blue
                  "Tumor"    = "#D55E00")   # orange

scRNA1$sample_type <- factor(scRNA1$sample_group, levels = c("Adjacent", "Tumor"))

p_umap <- DimPlot(
  scRNA1,
  group.by   = "sample_type",
  label      = FALSE,
  label.size = 5,
  reduction  = "umap",
  cols       = unname(pal_contrast[levels(scRNA1$sample_type)])
)
ggsave("sample_type_UMAP.pdf", plot = p_umap, width = 5, height = 4)
ggsave("sample_type_UMAP.png", plot = p_umap, width = 5, height = 4, dpi = 300, bg = "white")

message("✅ Fig 1B saved: sample_type_UMAP.pdf")

# ------------------------------------------------------------------
# 2. Fig 1D – Rain-cloud (box + density) for YBX1 expr by sample type
#    Helper function (reused by downstream scripts for other metrics)
# ------------------------------------------------------------------
plot_boxdensity_template <- function(df, out_prefix) {
  p.val <- kruskal.test(AUC ~ Subtype, data = df)
  p.lab <- paste0("P",
                  ifelse(p.val$p.value < 0.001, " < 0.001",
                         paste0(" = ", round(p.val$p.value, 3))))

  col1 <- "#0072B2"   # Adjacent
  col2 <- "#D55E00"   # Tumor
  col3 <- "#CC79A7"   # placeholder

  p_top <- ggplot(df, aes(x = AUC, color = Subtype, fill = Subtype)) +
    geom_density() +
    scale_color_manual(values = c(scales::alpha(col1, 0.7),
                                  scales::alpha(col2, 0.7),
                                  scales::alpha(col3, 0.7))) +
    scale_fill_manual(values  = c(scales::alpha(col1, 0.7),
                                  scales::alpha(col2, 0.7),
                                  scales::alpha(col3, 0.7))) +
    theme_classic() +
    xlab(paste0("Estimated AUC of ", unique(df$Drug))) +
    ylab(NULL) +
    theme(legend.position = "none",
          legend.title = element_blank(),
          axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y  = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_rug()

  p_bot <- ggplot(df, aes(Subtype, AUC, fill = Subtype)) +
    geom_boxplot(aes(col = Subtype)) +
    scale_fill_manual(values  = c(col1, col2, col3)) +
    scale_color_manual(values = c(col1, col2, col3)) +
    xlab(NULL) + ylab("Estimated AUC") +
    theme_void() +
    theme(legend.position = "right",
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 11, color = "black"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    annotate(geom = "text",
             x = 1.5, hjust = 1,
             y = max(df$AUC),
             size = 4, angle = 270, fontface = "bold",
             label = p.lab) +
    coord_flip()

  dat <- ggplot_build(p_bot)$data[[1]]
  p_bot <- p_bot + geom_segment(
    data = dat, aes(x = xmin, xend = xmax, y = middle, yend = middle),
    color = "white", inherit.aes = FALSE
  )

  p_all <- p_top %>% insert_bottom(p_bot, height = 0.4)
  ggsave(paste0(out_prefix, ".pdf"), p_all, width = 6, height = 3)
  ggsave(paste0(out_prefix, ".png"), p_all, width = 6, height = 3, dpi = 300, bg = "white")
  invisible(p_all)
}

message("✅ Fig 1D helper defined. Call plot_boxdensity_template() with YBX1 data to generate Fig 1D.")

# ------------------------------------------------------------------
# 3. Fig S1C – Stacked cell-type proportion bars
# ------------------------------------------------------------------
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

cell_order <- c(
  "B cell", "CD8+ T cell", "Cholangiocytes", "Dendritic cells", "Endothelial cell",
  "Fibroblast", "Hepatocytes", "moDCs", "Mono/Macro", "NK cells",
  "Plasma cell", "Plasmacytoid dendritic cells", "T cells"
)

need_cols <- c("curate_v1", "orig.ident")
miss <- setdiff(need_cols, colnames(scRNA1@meta.data))
if (length(miss)) stop("meta.data missing columns: ", paste(miss, collapse = ", "))

scRNA1$curate_v1 <- factor(as.character(scRNA1$curate_v1), levels = cell_order)
cells_keep  <- rownames(scRNA1@meta.data)[!is.na(scRNA1$curate_v1)]
scRNA1_use  <- subset(scRNA1, cells = cells_keep)

sample_order_given <- c(
  "GSM4116579_ICC_18_Adjacent", "GSM4116580_ICC_18_Tumor",
  "GSM4116581_ICC_20_Tumor",    "GSM4116582_ICC_23_Adjacent",
  "GSM4116583_ICC_23_Tumor",    "GSM4116584_ICC_24_Tumor1",
  "GSM4116585_ICC_24_Tumor2",   "GSM4116586_ICC_25_Adjacent"
)
present_samples <- intersect(sample_order_given, unique(scRNA1_use$orig.ident))
if (length(present_samples) == 0) present_samples <- unique(scRNA1_use$orig.ident)
scRNA1_use$orig.ident <- factor(as.character(scRNA1_use$orig.ident), levels = present_samples)

meta_df <- scRNA1_use@meta.data %>%
  dplyr::select(orig.ident, curate_v1)

count_df <- meta_df %>%
  count(orig.ident, curate_v1, name = "n_cells") %>%
  complete(orig.ident, curate_v1, fill = list(n_cells = 0)) %>%
  group_by(orig.ident) %>%
  mutate(total_cells = sum(n_cells),
         percent = 100 * n_cells / total_cells) %>%
  ungroup()

write.csv(count_df, "CellType_Proportion_by_Sample.csv", row.names = FALSE)

cell_present <- levels(droplevels(scRNA1_use$curate_v1))
pal_base     <- zzm60colors[seq_len(length(cell_present))]
names(pal_base) <- cell_present
pal_use <- colorspace::lighten(pal_base, amount = 0.18)
names(pal_use) <- names(pal_base)

# Per-sample stacked bar
p_stack <- ggplot(count_df, aes(x = orig.ident, y = percent, fill = curate_v1)) +
  geom_col(width = 0.9, color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.02))) +
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

ggsave("Stacked_Proportion_by_Sample.pdf", p_stack, width = 10, height = 4.8)
ggsave("Stacked_Proportion_by_Sample.png", p_stack, width = 10, height = 4.8, dpi = 300)

# Adjacent vs Tumor merged (Fig S1C)
AT_from_sample <- function(x) {
  ifelse(grepl("Adjacent|Normal|Para", x, ignore.case = TRUE), "Adjacent", "Tumor")
}

count_AT <- count_df %>%
  mutate(AT = factor(AT_from_sample(as.character(orig.ident)), levels = c("Adjacent", "Tumor"))) %>%
  group_by(AT, curate_v1) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop_last") %>%
  group_by(AT) %>%
  mutate(total_cells = sum(n_cells),
         percent = 100 * n_cells / total_cells) %>%
  ungroup()

p_stack_AT <- ggplot(count_AT, aes(x = AT, y = percent, fill = curate_v1)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = pal_use, name = "Cell type") +
  labs(x = NULL, y = "Proportion (%)",
       title = "Cell-type composition: Adjacent vs Tumor") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10)
  )

ggsave("Stacked_Proportion_AT.pdf",  p_stack_AT, width = 5.6, height = 4.4)
ggsave("Stacked_Proportion_AT.png",  p_stack_AT, width = 5.6, height = 4.4, dpi = 300)

message("✅ Fig S1C saved: Stacked_Proportion_AT.pdf")
message("✅ Per-sample stacked bars saved: Stacked_Proportion_by_Sample.pdf")
