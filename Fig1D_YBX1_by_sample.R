## ============================================================
## Fig 1D — YB-1 by sample origin (Tumor vs Adjacent), raincloud
## ------------------------------------------------------------
## Input : scRNA1   (loaded by 00_load_data.R; meta col sample_group, gene YBX1)
## Output: results/Fig1D_BoxDensity_YBX1_by_sample.pdf / .png
## Style : identical layout to Fig 1G (title "YB-1", coloured-text legend
##         top-right, two spread box plots, narrow). Colours per paper:
##         Adjacent=blue, Tumor=orange. Stats/data unchanged.
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(aplot)
  library(scales)
})

if (!exists("RESULTS")) {
  RESULTS <- Sys.getenv("YB1_RESULTS", unset = "../results")
  dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
}

stopifnot(exists("scRNA1"), inherits(scRNA1, "Seurat"))

# 与文章原图一致：Adjacent=蓝, Tumor=橙
pal_light <- c("Adjacent" = "#5AAFD6", "Tumor" = "#F29C64")
lev <- c("Adjacent", "Tumor")

vals <- FetchData(scRNA1, vars = "YBX1")[, 1]
grp  <- factor(as.character(scRNA1$sample_group), levels = lev)
df <- data.frame(Subtype = grp, AUC = as.numeric(vals))
df <- df[is.finite(df$AUC) & !is.na(df$Subtype), , drop = FALSE]

p.val <- tryCatch(wilcox.test(AUC ~ Subtype, data = df)$p.value, error = function(e) NA_real_)
p.lab <- if (is.finite(p.val)) paste0("P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ", round(p.val, 3)))) else "P = NA"

fill_vals  <- pal_light[lev]
alpha_vals <- alpha(fill_vals, 0.7); names(alpha_vals) <- lev

# 顶部密度：标题 "YB-1"，右上角彩色文字图例（Tumor 橙在上、Adjacent 蓝在下）
p_top <- ggplot(df, aes(x = AUC, color = Subtype, fill = Subtype)) +
  geom_density() +
  scale_color_manual(values = alpha_vals, na.translate = FALSE) +
  scale_fill_manual(values  = alpha_vals, na.translate = FALSE) +
  labs(title = "YB-1", x = NULL, y = NULL) +
  annotate("text", x = Inf, y = Inf, label = "Tumor",    color = pal_light[["Tumor"]],
           hjust = 1.1, vjust = 1.6, size = 4.6) +
  annotate("text", x = Inf, y = Inf, label = "Adjacent", color = pal_light[["Adjacent"]],
           hjust = 1.1, vjust = 3.0, size = 4.6) +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(4, 30, 2, 4)) +
  geom_rug()

# 底部箱线：两条岔开（Tumor 橙在上、Adjacent 蓝在下），无 y 标签
p_bot <- ggplot(df, aes(Subtype, AUC, fill = Subtype)) +
  geom_boxplot(aes(col = Subtype)) +
  scale_fill_manual(values = fill_vals, na.translate = FALSE) +
  scale_color_manual(values = fill_vals, na.translate = FALSE) +
  xlab(NULL) + ylab(NULL) + theme_void() +
  theme(legend.position = "none", legend.title = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = margin(2, 30, 2, 4)) +
  annotate(geom = "text", x = 1.5, hjust = 0.5, y = max(df$AUC, na.rm = TRUE) * 0.97,
           size = 3.3, angle = 270, fontface = "bold", label = p.lab) +
  coord_flip(clip = "off")

dat <- ggplot_build(p_bot)$data[[1]]
if (!is.null(dat) && nrow(dat)) {
  p_bot <- p_bot + geom_segment(data = dat, aes(x = xmin, xend = xmax, y = middle, yend = middle),
                                color = "white", inherit.aes = FALSE)
}

p_all <- p_top %>% insert_bottom(p_bot, height = 0.45)
ggsave(file.path(RESULTS, "Fig1D_BoxDensity_YBX1_by_sample.pdf"), p_all, width = 3.2, height = 3.4)
ggsave(file.path(RESULTS, "Fig1D_BoxDensity_YBX1_by_sample.png"), p_all, width = 3.2, height = 3.4, dpi = 300, bg = "white")
