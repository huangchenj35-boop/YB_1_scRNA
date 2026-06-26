## ============================================================
## Fig 1C — YB-1 density UMAP (Nebulosa weighted KDE)
## ------------------------------------------------------------
## Input : scRNA1   (loaded by 00_load_data.R) — FULL object (all cells),
##         gene YBX1, reduction umap. Fig 1C is on the total landscape
##         (not the epithelial subset used by Fig 1E / Fig 1F).
## Output: results/Fig1C_YBX1_density_UMAP.pdf / .png
## Match : paper uses Nebulosa::plot_density (viridis, legend "Density",
##         title "YB-1"). Primary path calls Nebulosa::plot_density.
##         Fallback reproduces Nebulosa's "ks" weighted-KDE directly when
##         the installed Seurat is too new for Nebulosa's FetchData() call
##         (defunct `slot=` arg in SeuratObject >= 5.0); identical math.
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Nebulosa)
  library(ggplot2)
})

if (!exists("RESULTS")) {
  RESULTS <- Sys.getenv("YB1_RESULTS", unset = "../results")
  dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
}

stopifnot(exists("scRNA1"), inherits(scRNA1, "Seurat"))

p <- tryCatch(
  plot_density(scRNA1, "YBX1", reduction = "umap") + ggtitle("YB-1"),
  error = function(e) {
    message("Nebulosa::plot_density unavailable on this Seurat version; ",
            "using its 'ks' weighted-KDE directly. (", conditionMessage(e), ")")
    suppressPackageStartupMessages(library(ks))
    emb <- Embeddings(scRNA1, "umap")[, 1:2]
    w   <- as.numeric(GetAssayData(scRNA1, assay = "RNA", layer = "data")["YBX1", ])
    wn  <- w / sum(w) * length(w)                       # Nebulosa weight normalisation
    H   <- ks::Hpi(emb, binned = TRUE)                  # Nebulosa "ks" bandwidth
    dens <- ks::kde(emb, w = wn, H = H, eval.points = emb)$estimate
    df <- data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2], Density = dens)
    df <- df[order(df$Density), ]                       # high density on top
    ggplot(df, aes(UMAP_1, UMAP_2, color = Density)) +
      geom_point(size = 0.4) +
      scale_color_viridis_c(name = "Density") +
      ggtitle("YB-1") + theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))
  }
)

ggsave(file.path(RESULTS, "Fig1C_YBX1_density_UMAP.pdf"), p, width = 7, height = 6)
ggsave(file.path(RESULTS, "Fig1C_YBX1_density_UMAP.png"), p, width = 7, height = 6, dpi = 300, bg = "white")
