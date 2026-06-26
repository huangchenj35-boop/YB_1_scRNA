## ============================================================
## Fig S1F — Cholangiocyte subtype marker stacked violin
## ------------------------------------------------------------
## Source : 03_Script82_cholangiocyte_subclustering.R
##          - cluster -> subtype mapping (lines 13-45)
##          - long-format violin table vln.df (lines 47-83)
##          The original file stops at `head(vln.df)`; the violin ggplot
##          (truncated in the supplied script) is completed here to match
##          the article Fig S1F (10 markers x 3 subtypes, horizontal stack).
## Input  : scRNA_sub (loaded by 00_load_data.R; seurat_clusters, RNA assay)
## Output : results/FigS1F_subtype_marker_violin.pdf / .png
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)    # 原脚本用 library(tidyverse)，此处仅取其中的 dplyr + tidyr
  library(tidyr)
  library(ggplot2)
})

if (!exists("RESULTS")) {
  RESULTS <- Sys.getenv("YB1_RESULTS", unset = "../results")
  dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
}

stopifnot(exists("scRNA_sub"), "seurat_clusters" %in% colnames(scRNA_sub@meta.data))
obj <- scRNA_sub
DefaultAssay(obj) <- "RNA"

## ---- 1) cluster -> subtype (verbatim from 03_Script82), paper labels ----
cluster2type <- c(
  `0`="Malignant", `1`="Malignant", `2`="Cholangiocyte", `3`="Malignant",
  `4`="Malignant", `5`="Cholangiocyte", `6`="Cholangiocyte", `7`="Cholangiocyte",
  `8`="Cholangiocyte", `9`="Cholangiocyte", `10`="Malignant", `11`="LPC",
  `12`="Malignant", `13`="Malignant"
)
clusters <- as.character(obj@meta.data$seurat_clusters)
obj@meta.data$curate_v1 <- factor(unname(cluster2type[clusters]),
                                  levels = c("LPC", "Cholangiocyte", "Malignant"))

## ---- 2) long table vln.df (as in 03_Script82), article's 10 markers ----
features <- c("HNF4A","ALB","SOX9","KRT19","MYC","BIRC5","UBE2C","CDK1","CCNB1","CCNB2")
features <- features[features %in% rownames(obj[["RNA"]])]

df <- FetchData(obj, vars = c(features, "curate_v1"))
df$celltype <- df$curate_v1
vln.df <- df %>%
  tidyr::pivot_longer(cols = all_of(features), names_to = "gene", values_to = "exp") %>%
  dplyr::mutate(gene = factor(gene, levels = features),
                celltype = factor(celltype, levels = c("LPC","Cholangiocyte","Malignant")))

## ---- 3) completed stacked violin matching article Fig S1F ----
cols <- c("LPC" = "#185070", "Cholangiocyte" = "#809070", "Malignant" = "#F0D868")

p <- ggplot(vln.df, aes(x = celltype, y = exp, fill = celltype)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
  facet_grid(rows = vars(gene)) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.spacing = unit(0, "lines"),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x.top = element_text(angle = 45, hjust = 0, color = "black"),
    strip.text.y = element_text(angle = 0, hjust = 0),
    strip.background = element_blank()
  )

ggsave(file.path(RESULTS, "FigS1F_subtype_marker_violin.pdf"), plot = p, width = 4.5, height = 6)
ggsave(file.path(RESULTS, "FigS1F_subtype_marker_violin.png"), plot = p, width = 4.5, height = 6, dpi = 300, bg = "white")
