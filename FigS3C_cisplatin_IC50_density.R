## ============================================================
## Fig S3C — Cisplatin lnIC50 density on UMAP (epithelial subset)
## ------------------------------------------------------------
## Source : 16_Script21_FigS3C_UMAP.R — ensure_meta_as_assay_feature() +
##          Nebulosa::plot_density(); ported verbatim. Input object is the
##          author's own scRNAauc_with_CisplatinGroup_Cisplatin.rds, which
##          already carries lnIC50_Cisplatin_1005 on the figure's UMAP.
## Input  : scRNAauc_cis  (loaded by 00_load_data.R)
## Output : results/FigS3C_cisplatin_lnIC50_density_UMAP.pdf / .png
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Nebulosa)
  library(Matrix)
  library(ggplot2)
})

if (!exists("RESULTS")) {
  RESULTS <- Sys.getenv("YB1_RESULTS", unset = "../results")
  dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
}

stopifnot(exists("scRNAauc_cis"), inherits(scRNAauc_cis, "Seurat"),
          "lnIC50_Cisplatin_1005" %in% colnames(scRNAauc_cis@meta.data))
ln_ic_col <- "lnIC50_Cisplatin_1005"
red_use   <- "umap"   # 文章 S3C 标题 "on UMAP"，用 umap（非 tsne）

## ---- 16_Script21 source: put a meta column into a META assay as a feature ----
ensure_meta_as_assay_feature <- function(obj, meta_col, assay_name = "META"){
  stopifnot(meta_col %in% colnames(obj@meta.data))
  v <- obj@meta.data[[meta_col]]; v[!is.finite(v)] <- 0
  cells <- colnames(obj)
  m_counts <- sparseMatrix(
    i = rep(1L, length(cells)),
    j = seq_along(cells),
    x = as.numeric(v),
    dims = c(1L, length(cells)),
    dimnames = list(meta_col, cells)
  )
  if (!(assay_name %in% names(obj@assays))){
    assay <- CreateAssayObject(counts = m_counts)
    assay@data <- m_counts
    obj[[assay_name]] <- assay
  } else {
    old <- GetAssayData(obj, assay = assay_name, slot = "counts")
    if (meta_col %in% rownames(old)) old[meta_col, ] <- m_counts[1, ] else old <- rbind(old, m_counts)
    obj@assays[[assay_name]]@counts <- old
    obj@assays[[assay_name]]@data   <- old
  }
  obj
}

ttl <- "lnIC50 Cisplatin 1005 density on UMAP"

p <- tryCatch({
  obj <- ensure_meta_as_assay_feature(scRNAauc_cis, ln_ic_col, assay_name = "META")
  DefaultAssay(obj) <- "META"
  Nebulosa::plot_density(object = obj, features = ln_ic_col, reduction = red_use) + ggtitle(ttl)
}, error = function(e) {
  message("Nebulosa::plot_density unavailable on this Seurat version; ",
          "using its 'ks' weighted-KDE directly. (", conditionMessage(e), ")")
  suppressPackageStartupMessages(library(ks))
  emb <- Embeddings(scRNAauc_cis, red_use)[, 1:2]
  w   <- scRNAauc_cis@meta.data[[ln_ic_col]]; w[!is.finite(w)] <- 0
  wn  <- w / sum(w) * length(w)
  H   <- ks::Hpi(emb, binned = TRUE)
  dens <- ks::kde(emb, w = wn, H = H, eval.points = emb)$estimate
  df <- data.frame(D1 = emb[, 1], D2 = emb[, 2], Density = dens)
  df <- df[order(df$Density), ]
  ggplot(df, aes(D1, D2, color = Density)) +
    geom_point(size = 0.4) +
    scale_color_viridis_c(name = "Density") +
    labs(title = ttl, x = "UMAP_1", y = "UMAP_2") +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5))
})

ggsave(file.path(RESULTS, "FigS3C_cisplatin_lnIC50_density_UMAP.pdf"), p, width = 7, height = 6)
ggsave(file.path(RESULTS, "FigS3C_cisplatin_lnIC50_density_UMAP.png"), p, width = 7, height = 6, dpi = 300, bg = "white")
