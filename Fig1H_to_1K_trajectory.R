## ============================================================
## Fig 1H-K — Monocle2 (DDRTree) pseudotime trajectory
## Source: 14_Script81_monocle2_pseudotime.R (ported verbatim).
## Adapted: input scRNA1.subset -> scRNA_sub; outputs -> results/ as Fig1H/I/J/K.
## Output: results/Fig1H_pseudotime.pdf, Fig1I_group_copykat.pdf,
##         Fig1J_state.pdf, Fig1K_YBX1_in_pseudotime.pdf
## ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(monocle)       # v2
  library(Matrix)
  library(Biobase)
  library(ggplot2)
  library(dplyr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

set.seed(1234)

## --- compatibility shim: monocle 2.x calls dplyr SE verbs (group_by_/summarize_/...)
## defunct in dplyr >= 1.1; restore them by delegating to the modern verbs.
suppressWarnings(suppressMessages({
  if (requireNamespace("rlang", quietly = TRUE)) {
    .se <- function(modern, named = FALSE) function(.data, ..., .dots) {
      d <- if (!missing(.dots)) .dots else list(...)
      qs <- lapply(d, function(x) if (is.character(x)) rlang::parse_expr(x)
                   else tryCatch(rlang::get_expr(rlang::as_quosure(x)), error = function(e) x))
      if (named) names(qs) <- names(d)
      do.call(modern, c(list(.data), qs))
    }
    ns <- asNamespace("dplyr")
    for (v in c("group_by_","arrange_","filter_","select_","distinct_"))
      try(assignInNamespace(v, .se(getExportedValue("dplyr", sub("_$","",v))), ns), silent = TRUE)
    for (v in c("summarize_","summarise_","mutate_"))
      try(assignInNamespace(v, .se(getExportedValue("dplyr", sub("_$","",v)), named = TRUE), ns), silent = TRUE)
  }
  ## monocle 2.x calls igraph::dfs(neimode=) — defunct in igraph >= 1.3; map to mode=
  if (requireNamespace("igraph", quietly = TRUE)) {
    ig <- asNamespace("igraph"); real_dfs <- get("dfs", ig)
    shim_dfs <- function(graph, root, neimode, mode, ...) {
      if (!missing(neimode) && missing(mode)) mode <- neimode
      a <- list(graph = graph, root = root, ...); if (!missing(mode)) a$mode <- mode
      do.call(real_dfs, a)
    }
    try(assignInNamespace("dfs", shim_dfs, ig), silent = TRUE)
  }
  ## monocle 2.x uses V(g)[nei(...)]; igraph >= 2.1 renamed nei() -> .nei().
  ## Patch the monocle functions that still call nei( to use .nei(.
  if (requireNamespace("monocle", quietly = TRUE) && requireNamespace("igraph", quietly = TRUE) &&
      exists(".nei", asNamespace("igraph"))) {
    mns <- asNamespace("monocle")
    for (f in ls(mns, all.names = TRUE)) {
      obj <- tryCatch(get(f, mns), error = function(e) NULL)
      if (!is.function(obj)) next
      src <- deparse(obj)
      if (any(grepl("[^.A-Za-z0-9_]nei\\(", src))) {
        src2 <- gsub("([^.A-Za-z0-9_])nei\\(", "\\1.nei(", src)
        newf <- tryCatch(eval(parse(text = paste(src2, collapse = "\n")), envir = mns),
                         error = function(e) NULL)
        if (is.function(newf)) { environment(newf) <- mns; try(assignInNamespace(f, newf, "monocle"), silent = TRUE) }
      }
    }
  }
}))

if (!exists("RESULTS")) { RESULTS <- Sys.getenv("YB1_RESULTS", unset="../results"); dir.create(RESULTS, recursive=TRUE, showWarnings=FALSE) }
stopifnot(exists("scRNA_sub"))
scRNAauc_Cisplatin <- scRNA_sub
## ========= 0) 前置检查 & 只取 diploid / aneuploid =========
stopifnot(exists("scRNAauc_Cisplatin"), inherits(scRNAauc_Cisplatin, "Seurat"))
if (!"group_copykat" %in% colnames(scRNAauc_Cisplatin@meta.data)) {
  stop("在 scRNAauc_Cisplatin@meta.data 中找不到列：group_copykat")
}
scRNAauc_Cisplatin$group_copykat <- factor(as.character(scRNAauc_Cisplatin$group_copykat),
                                           levels = c("diploid","aneuploid"))

.cdscache <- Sys.getenv("YB1_CDS_CACHE", "")
if (nzchar(.cdscache) && file.exists(.cdscache)) {
message("loaded cached cds from ", .cdscache); cds <- readRDS(.cdscache)
} else {
cells_use <- rownames(scRNAauc_Cisplatin@meta.data)[
  !is.na(scRNAauc_Cisplatin$group_copykat) &
    scRNAauc_Cisplatin$group_copykat %in% c("diploid","aneuploid")
]
if (length(cells_use) < 50) stop("可用细胞过少（<50）。")

obj <- subset(scRNAauc_Cisplatin, cells = cells_use)

## ========= 1) 取计数矩阵 & 构建 Monocle v2 CellDataSet =========
assay_use <- if ("RNA" %in% names(obj@assays)) "RNA" else DefaultAssay(obj)
cnt <- tryCatch(GetAssayData(obj, assay = assay_use, layer = "counts"), error=function(e) GetAssayData(obj, assay = assay_use, slot = "counts"))
if (is.null(cnt) || length(cnt) == 0) {
  message("⚠️ counts 槽为空，退回使用 data 槓。")
  cnt <- tryCatch(GetAssayData(obj, assay = assay_use, layer = "data"), error=function(e) GetAssayData(obj, assay = assay_use, slot = "data"))
}
if (!inherits(cnt, "dgCMatrix")) cnt <- as(as.matrix(cnt), "dgCMatrix")

# 构建 pheno/feature 数据
pd <- new("AnnotatedDataFrame", data = obj@meta.data[colnames(cnt), , drop = FALSE])

# gene_short_name：若行为 ENSEMBL，尝试映射为 SYMBOL；否则用行名
rn <- rownames(cnt)
if (any(grepl("^ENSG", rn))) {
  ann <- AnnotationDbi::select(org.Hs.eg.db,
                               keys    = rn,
                               keytype = "ENSEMBL",
                               columns = "SYMBOL")
  ann <- ann[!is.na(ann$SYMBOL), ]
  gene_short_name <- setNames(ann$SYMBOL, ann$ENSEMBL)[rn]
  gene_short_name[is.na(gene_short_name)] <- rn[is.na(gene_short_name)]
} else {
  gene_short_name <- rn
}
fd <- new("AnnotatedDataFrame",
          data = data.frame(gene_short_name = gene_short_name,
                            row.names = rn, check.names = FALSE))

cds <- newCellDataSet(
  cnt,
  phenoData            = pd,
  featureData          = fd,
  lowerDetectionLimit  = 0.5,                 # 保持你的参数
  expressionFamily     = negbinomial.size()   # 保持你的参数
)

# 尺度 & 离散
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)       # 保持你的参数
message("表达基因统计完成；cells=", ncol(cds), ", genes=", nrow(cds))

## ========= 2) 选择排序基因（ordering genes） =========
disp_table <- dispersionTable(cds)
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.1 &
                           dispersion_empirical >= 1 * dispersion_fit)$gene_id
if (length(ordering_genes) < 500) {
  # 兜底：加入 diploid vs aneuploid 的 DE 基因（简单检验），增加轨迹信息量
  message("⚠️ ordering genes 偏少，增加两组差异驱动基因作为补充。")
  grp <- factor(pData(cds)$group_copykat, levels = c("diploid","aneuploid"))
  mu1 <- Matrix::rowMeans(exprs(cds)[, grp=="diploid", drop=FALSE] > 0)
  mu2 <- Matrix::rowMeans(exprs(cds)[, grp=="aneuploid", drop=FALSE] > 0)
  de_hint <- names(sort(abs(mu2 - mu1), decreasing = TRUE))[1:2000]
  de_hint <- intersect(de_hint, rownames(cds))
  ordering_genes <- unique(c(ordering_genes, de_hint))
}
cds <- setOrderingFilter(cds, ordering_genes)

## ========= 3) 降维 + 拟时序排序（DDRTree） =========
cds <- reduceDimension(cds,
                       max_components = 2,
                       method = "DDRTree")    # 保持你的参数

## ========= 4) 设置二倍体细胞作为发育起点 =========
# 首先运行 orderCells() 不带 root_state 参数
cds <- orderCells(cds)

# 然后指定 diploid 作为发育起点
cds <- orderCells(cds, root_state = "5")
if (nzchar(.cdscache)) try(saveRDS(cds, .cdscache), silent = TRUE)
}  ## end build-or-load-cds

## ========= Custom trajectory plotting =========
## Draw the trajectory (Component 1/2 + DDRTree backbone) directly from the cds,
## bypassing monocle::plot_cell_trajectory (which calls dplyr SE verbs that are
## defunct on this R). Output is the same Monocle2 DDRTree trajectory.
suppressPackageStartupMessages({ library(ggplot2); library(igraph) })
S <- t(monocle::reducedDimS(cds)); colnames(S) <- c("C1", "C2")
K <- t(monocle::reducedDimK(cds)); rownames(K) <- colnames(monocle::reducedDimK(cds))
el  <- igraph::as_edgelist(monocle::minSpanningTree(cds))
seg <- data.frame(x = K[el[,1],1], y = K[el[,1],2], xend = K[el[,2],1], yend = K[el[,2],2])
## vertical flip so the trunk points up like the article (DDRTree orientation is arbitrary)
S[,2] <- -S[,2]; K[,2] <- -K[,2]
seg <- data.frame(x = K[el[,1],1], y = K[el[,1],2], xend = K[el[,2],1], yend = K[el[,2],2])
.gid <- intersect(c("YBX1","ybx1","Ybx1"), rownames(cds))[1]
## YB-1 colouring uses the Seurat normalised expression (0–6), like the article,
## not the raw counts stored in the monocle cds.
.ybnorm <- tryCatch(as.numeric(GetAssayData(scRNA_sub, assay = "RNA", layer = "data")[.gid, colnames(cds)]),
                    error = function(e) as.numeric(log1p(exprs(cds[.gid, , drop = FALSE]))))
pdat <- data.frame(C1 = S[,1], C2 = S[,2],
                   Pseudotime    = pData(cds)$Pseudotime,
                   State         = factor(pData(cds)$State),
                   group_copykat = pData(cds)$group_copykat,
                   YBX1          = .ybnorm)

traj_plot <- function(color_by, title, scale_layer) {
  ggplot() +
    geom_segment(data = seg, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.6, color = "black") +
    geom_point(data = pdat, aes(C1, C2, color = .data[[color_by]]), size = 0.6) +
    labs(title = title, x = "Component 1", y = "Component 2") +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + scale_layer
}
save2 <- function(p, base, w = 5.5, h = 4.5) {
  ggsave(file.path(RESULTS, paste0(base, ".pdf")), p, width = w, height = h)
  ggsave(file.path(RESULTS, paste0(base, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
}
## 与原图一致的配色
copy_pal  <- c("diploid" = "#F8766D", "aneuploid" = "#00BFC4")             # 原图: diploid 肉粉, aneuploid 青
state_pal <- c("1"="#C77CFF","2"="#FF61CC","3"="#A3A500","4"="#00BF7D","5"="#00B0F6")

save2(traj_plot("Pseudotime", "Pseudotime",
                scale_color_gradient(low = "#08306b", high = "#c6dbef", name = "Pseudotime")), "Fig1H_pseudotime")
save2(traj_plot("group_copykat", "CopyKAT",
                scale_color_manual(values = copy_pal, name = "")), "Fig1I_group_copykat")
save2(traj_plot("State", "State",
                scale_color_manual(values = state_pal, name = "State")), "Fig1J_state")
save2(traj_plot("YBX1", "YB-1",
                scale_color_gradientn(colours = c("#2166ac","#f7fbff","#b2182b"), name = "YB-1")), "Fig1K_YBX1")
message("DONE: Fig1H/1I/1J/1K trajectory plots written to ", RESULTS)
