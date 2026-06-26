## ============================================================
## Fig 1G — YB-1 box-density by CopyKAT group (aneuploid vs diploid)
## ------------------------------------------------------------
## Input : scRNA_sub (loaded by 00_load_data.R) — epithelial subset, meta col
##         group_copykat (aneuploid/diploid), gene YBX1. (Fig 1D uses full scRNA1.)
## Output: results/Fig1G_BoxDensity_YBX1_by_group_copykat_*.pdf / .png
## Note  : This is also the shared box-density engine for Fig 3B/3C (when SCENIC
##         AUC/GSVA columns are present) and Fig 3D/3E (when run on scRNA_cis /
##         scRNA_gem). Plotting logic is unchanged; only header/guard added.
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)    # 仍保留，接口不改
  library(aplot)        # 用 insert_bottom 组合密度+箱线
  library(ggpubr)       # 计算/标注 p 值（备用）
  library(scales)       # alpha()
})

# allow standalone use if 00_load_data.R was not sourced
if (!exists("RESULTS")) {
  RESULTS <- Sys.getenv("YB1_RESULTS", unset = "../results")
  dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
}

`%||%` <- function(a,b) if (!is.null(a) && length(a)>0) a else b
safe <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

## ========= 前置：确保对象 =========
## 用 scRNAauc(= 上皮子集 12290，已含 YBX1_extended_907g 与 GSVA 列）：
## 同一套 raincloud 引擎一次产出 Fig1G(YBX1)+ Fig3B下(regulon AUC)+ Fig3C下(GSVA)。
## 若 scRNAauc 不在(未提供 SCENIC 文件)，回退用 scRNA_sub 只出 Fig1G。
obj_src <- if (exists("scRNAauc") && inherits(scRNAauc, "Seurat")) scRNAauc else scRNA_sub
stopifnot(inherits(obj_src, "Seurat"), "group_copykat" %in% colnames(obj_src@meta.data))

## ========= 小工具 =========
# 找 YBX1 基因（仅用“原始表达”画图，不对表达做z）
find_ybx1_gene <- function(obj){
  assays <- names(obj@assays)
  assays <- unique(c(intersect(c("RNA","SCT"), assays), setdiff(assays, c("RNA","SCT"))))
  for (a in assays) {
    rn <- tryCatch(rownames(obj@assays[[a]]), error = function(e) character(0))
    if (length(rn) == 0) next
    if ("YBX1" %in% rn) return(list(assay=a, feature="YBX1"))
    hit <- grep("^YBX1$", rn, ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(list(assay=a, feature=hit[1]))
  }
  stop("未在任何 assay 的行名中找到 YBX1。")
}

# 取/写入 YBX1 AUC（若当前对象没有则尝试从 scRNAauc 抄）
get_ybx1_auc <- function(obj){
  mdn <- colnames(obj@meta.data)
  cand <- grep("^YBX1", mdn, ignore.case = TRUE, value = TRUE)
  cand <- cand[!cand %in% c("YBX1","YBX1_expr_z")]  # 排除表达/旧 z
  if (length(cand)) {
    non_ext <- cand[!grepl("_extended$", cand, ignore.case = TRUE)]
    return(list(obj=obj, col = non_ext[1] %||% cand[1]))
  }
  if (exists("scRNAauc") && inherits(scRNAauc, "Seurat")) {
    hit <- do.call(rbind, lapply(names(scRNAauc@assays), function(a){
      rn <- tryCatch(rownames(scRNAauc@assays[[a]]), error = function(e) character(0))
      if (!length(rn)) return(NULL)
      idx <- grep("^YBX1", rn, ignore.case = TRUE, value = TRUE)
      if (length(idx)) data.frame(assay=a, feature=idx, stringsAsFactors=FALSE) else NULL
    }))
    if (!is.null(hit) && nrow(hit)) {
      non_ext <- hit[!grepl("_extended$", hit$feature, ignore.case = TRUE), , drop = FALSE]
      pick <- if (nrow(non_ext)) non_ext[1, ] else hit[1, ]
      auc_vec <- as.numeric(GetAssayData(scRNAauc, assay = pick$assay, slot = "data")[pick$feature, ])
      names(auc_vec) <- colnames(scRNAauc)
      v <- rep(NA_real_, ncol(obj)); names(v) <- colnames(obj)
      common <- intersect(names(auc_vec), names(v))
      v[common] <- auc_vec[common]
      new_col <- gsub("[ ()]", "_", pick$feature) # 安全列名
      obj[[new_col]] <- v
      return(list(obj=obj, col=new_col))
    }
  }
  list(obj=obj, col=NA_character_)
}

## —— 确保 YBX1 targets GSVA 列（与你前面写入 meta 的列名一致），
## —— 在子集对象缺列时按细胞名从全量 scRNA1 对齐拷贝
ensure_targets_cols <- function(obj){
  target_raw <- "GSVA_YBX1_targets_ssgsea"
  target_z   <- "GSVA_YBX1_targets_ssgsea_z"
  
  copy_if_missing <- function(dst, src, col) {
    if (!(col %in% colnames(dst@meta.data)) && (col %in% colnames(src@meta.data))) {
      v <- rep(NA_real_, ncol(dst)); names(v) <- colnames(dst)
      cc <- intersect(names(v), rownames(src@meta.data))
      v[cc] <- src@meta.data[cc, col]
      dst[[col]] <- v
    }
    dst
  }
  if (exists("scRNA1") && inherits(scRNA1, "Seurat")) {
    obj <- copy_if_missing(obj, scRNA1, target_raw)
    obj <- copy_if_missing(obj, scRNA1, target_z)
  }
  
  mdn <- colnames(obj@meta.data)
  raw_col <- if (target_raw %in% mdn) target_raw else NA_character_
  z_col   <- if (target_z   %in% mdn) target_z   else NA_character_
  
  # 只有 raw 没有 z：现算 z
  if (!is.na(raw_col) && is.na(z_col)) {
    vv <- obj@meta.data[[raw_col]]
    mu <- mean(vv, na.rm=TRUE); sdv <- stats::sd(vv, na.rm=TRUE)
    z_col <- paste0(raw_col, "_z")
    obj[[z_col]] <- if (is.finite(sdv) && sdv>0) (vv-mu)/sdv else rep(0, length(vv))
  }
  # 只有 z 没有 raw：raw 用 z 代画
  if (is.na(raw_col) && !is.na(z_col)) raw_col <- z_col
  
  # 兜底模糊匹配
  if (is.na(raw_col) && is.na(z_col)) {
    cand_raw <- grep("^GSVA[_-]?YBX1.*targets.*ssgsea$", mdn, ignore.case = TRUE, value = TRUE)
    cand_z   <- grep("^GSVA[_-]?YBX1.*targets.*ssgsea[_-]?z$", mdn, ignore.case = TRUE, value = TRUE)
    raw_col <- if (length(cand_raw)) cand_raw[1] else NA_character_
    z_col   <- if (length(cand_z))   cand_z[1]   else NA_character_
    if (!is.na(raw_col) && is.na(z_col)) {
      vv <- obj@meta.data[[raw_col]]
      mu <- mean(vv, na.rm=TRUE); sdv <- stats::sd(vv, na.rm=TRUE)
      z_col <- paste0(raw_col, "_z")
      obj[[z_col]] <- if (is.finite(sdv) && sdv>0) (vv-mu)/sdv else rep(0, length(vv))
    }
    if (is.na(raw_col) && !is.na(z_col)) raw_col <- z_col
  }
  list(obj=obj, raw=raw_col, z=z_col)
}

# 深色同色系配色（只改这里）
pal_light <- c(
  "aneuploid" = "#F29C64",  # 橙 — 与文章原图一致 (paper Fig1G: aneuploid=orange)
  "diploid"   = "#5AAFD6"   # 蓝 — 与文章原图一致 (paper Fig1G: diploid=blue)
)

## ========= 峰峦替换为：密度 + 箱线（FigureYa227 风格）；其余参数不变 =========
plot_and_save <- function(obj, feature, prefix, group_col, tag){
  # 准备数据（分组=group_copykat）
  stopifnot(group_col %in% colnames(obj@meta.data))
  vals <- FetchData(obj, vars = feature)[,1]
  grp  <- obj@meta.data[[group_col]]
  # 固定顺序：aneuploid -> diploid
  grp  <- factor(as.character(grp), levels = c("aneuploid","diploid"))
  df <- data.frame(
    Subtype = grp,
    AUC     = as.numeric(vals),
    Drug    = prefix,
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$AUC) & !is.na(df$Subtype), , drop = FALSE]
  if (nrow(df) == 0) stop("没有有效数值用于绘图：", prefix)
  
  # 计算 P 值：两组用 wilcox，>2 组 kruskal；只有一组则 NA
  if (nlevels(df$Subtype) < 2) {
    p.val <- NA_real_
  } else if (nlevels(df$Subtype) == 2) {
    p.val <- tryCatch(wilcox.test(AUC ~ Subtype, data = df)$p.value, error = function(e) NA_real_)
  } else {
    p.val <- tryCatch(kruskal.test(AUC ~ Subtype, data = df)$p.value, error = function(e) NA_real_)
  }
  p.lab <- if (is.finite(p.val)) {
    paste0("P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ", round(p.val, 3))))
  } else "P = NA"
  
  # 颜色（顶层透明 0.7）
  fill_vals  <- pal_light[levels(df$Subtype)]
  alpha_vals <- alpha(fill_vals, 0.7)
  names(alpha_vals) <- levels(df$Subtype)
  
  # 顶部密度（标题居中；图例改为右上角彩色文字，与文章原图一致）
  title_lab <- if (toupper(prefix) == "YBX1") "YB-1" else prefix
  lev <- levels(df$Subtype)
  p_top <- ggplot(df, aes(x = AUC, color = Subtype, fill = Subtype)) +
    geom_density() +
    scale_color_manual(values = alpha_vals, na.translate = FALSE) +
    scale_fill_manual(values  = alpha_vals, na.translate = FALSE) +
    labs(title = title_lab, x = NULL, y = NULL) +
    annotate("text", x = Inf, y = Inf, label = lev[1], color = fill_vals[[lev[1]]],
             hjust = 1.1, vjust = 1.6, size = 4.6) +
    annotate("text", x = Inf, y = Inf, label = lev[2], color = fill_vals[[lev[2]]],
             hjust = 1.1, vjust = 3.0, size = 4.6) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 30, 2, 4)
    ) +
    geom_rug()
  
  # 底部箱线
  p_bot <- ggplot(df, aes(Subtype, AUC, fill = Subtype)) +
    geom_boxplot(aes(col = Subtype)) +
    scale_fill_manual(values = fill_vals, na.translate = FALSE) +
    scale_color_manual(values = fill_vals, na.translate = FALSE) +
    xlab(NULL) + ylab(NULL) +
    theme_void() +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(2, 30, 2, 4)
    ) +
    scale_x_discrete(limits = rev(lev)) +   # 文章原图：aneuploid 在上、diploid 在下
    annotate(
      geom="text", x = 1.5, hjust = 0.5, y = max(df$AUC, na.rm = TRUE) * 0.97,
      size = 3.3, angle = 270, fontface = "bold", label = p.lab
    ) +
    coord_flip(clip = "off")
  
  # 用白色标记箱子的中位数
  dat <- ggplot_build(p_bot)$data[[1]]
  if (!is.null(dat) && nrow(dat)) {
    p_bot <- p_bot + geom_segment(
      data = dat,
      aes(x = xmin, xend = xmax, y = middle, yend = middle),
      color = "white", inherit.aes = FALSE
    )
  }
  
  # 组合 & 保存
  p_all <- p_top %>% insert_bottom(p_bot, height = 0.45)
  
  # output naming only: map this panel to its article figure id
  fig_tag <- if (grepl("Cisplatin",  group_col, ignore.case = TRUE)) "Fig3D_" else
             if (grepl("Gemcitabine", group_col, ignore.case = TRUE)) "Fig3E_" else
             if (identical(prefix, "YBX1"))                            "Fig1G_" else
             if (grepl("907|AUC|regulon|extended", prefix, ignore.case = TRUE)) "Fig3Bbottom_" else
             if (grepl("GSVA|target", prefix, ignore.case = TRUE))     "Fig3Cbottom_" else ""
  fn_base <- file.path(RESULTS,
                       paste0(fig_tag, "BoxDensity_", safe(prefix), "_by_", group_col, "_", tag))
  ggsave(paste0(fn_base, ".pdf"), p_all, width = 3.2, height = 3.4)
  ggsave(paste0(fn_base, ".png"), p_all, width = 3.2, height = 3.4, dpi = 300, bg = "white")
  
  invisible(list(p_top = p_top, p_bot = p_bot, p_all = p_all))
}

## ========= 主流程（按二倍体/非二倍体分组画“三套图”） =========
run_set <- function(obj, tag){
  # 固定分组列：group_copykat（aneuploid / diploid）
  obj$group_copykat <- factor(as.character(obj$group_copykat), levels = c("aneuploid","diploid"))
  
  # 找 YBX1（只画原始表达）
  yinfo <- find_ybx1_gene(obj)
  DefaultAssay(obj) <- yinfo$assay
  feat_y <- yinfo$feature
  
  # YBX1 AUC（raw+z；可从 scRNAauc 拿）
  auc_info <- get_ybx1_auc(obj); obj <- auc_info$obj; auc_col <- auc_info$col
  if (!is.na(auc_col) && nzchar(auc_col)) {
    z_auc <- paste0(auc_col, "_z")
    if (!(z_auc %in% colnames(obj@meta.data))) {
      vv <- obj@meta.data[[auc_col]]
      mu <- mean(vv, na.rm=TRUE); sdv <- stats::sd(vv, na.rm=TRUE)
      obj[[z_auc]] <- if (is.finite(sdv) && sdv>0) (vv-mu)/sdv else rep(0, length(vv))
    }
  }
  
  # YBX1 targets（raw+z）
  tg_info <- ensure_targets_cols(obj); obj <- tg_info$obj
  gsva_raw_col <- tg_info$raw; gsva_z_col <- tg_info$z
  
  grp <- "group_copykat"  # ← 分组列
  
  # 1) YBX1 原始表达
  plot_and_save(obj, feat_y, "YBX1", grp, tag)
  
  # 2) YBX1 AUC（若可用，raw+z）
  if (!is.na(auc_col) && nzchar(auc_col)) {
    plot_and_save(obj, auc_col,               safe(auc_col),        grp, tag)
    plot_and_save(obj, paste0(auc_col, "_z"), safe(paste0(auc_col, "_z")), grp, tag)
  } else {
    message("⚠️ [", tag, "] 未找到 YBX1 AUC 列，跳过 AUC 绘图。")
  }
  
  # 3) Targets GSVA（若可用，raw+z）
  if (!is.na(gsva_raw_col)) plot_and_save(obj, gsva_raw_col, safe(gsva_raw_col), grp, tag) else
    message("⚠️ [", tag, "] 未找到 YBX1 targets GSVA（raw）列。")
  if (!is.na(gsva_z_col))   plot_and_save(obj, gsva_z_col,   safe(gsva_z_col),   grp, tag) else
    message("⚠️ [", tag, "] 未找到 YBX1 targets GSVA（z）列。")
  
  invisible(obj)
}

## —— 仅画“全体细胞”（包含 aneuploid & diploid 两组）—— ##
obj_all <- obj_src     # scRNAauc(优先) 或 scRNA_sub:上皮子集 aneuploid 7040 / diploid 5250
obj_all <- run_set(obj_all, tag = "AllCells_group_copykat")

message("✅ 完成：已按 aneuploid vs diploid 输出三套图（YBX1、YBX1 AUC、YBX1 targets；AUC/targets 各 raw+z）到 output/，配色为更深同色系。")

