## ============================================================
## Fig 3F — ABCB1 / ABCC1 / ABCC2 / MVP density UMAP (full object)
## ------------------------------------------------------------
## Source : Nebulosa::plot_density density UMAP — same plotting routine the
##          authors use in 16_Script21_FigS3C_UMAP.R, here applied to the four
##          ABC-transporter / vault genes. Only the feature names and (full)
##          object differ; the density routine is unchanged.
## Input  : scRNA1 (full object, reduction umap; genes ABCB1/ABCC1/ABCC2/MVP)
## Output : results/Fig3F_<gene>_density_UMAP.png  +  combined row
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Nebulosa)
  library(ggplot2)
  library(patchwork)
})

if (!exists("RESULTS")) {
  RESULTS <- Sys.getenv("YB1_RESULTS", unset = "../results")
  dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
}

stopifnot(exists("scRNA1"), inherits(scRNA1, "Seurat"))

genes <- c("ABCB1", "ABCC1", "ABCC2", "MVP")
genes <- genes[genes %in% rownames(scRNA1[["RNA"]])]

density_one <- function(obj, g) {
  tryCatch(
    plot_density(obj, g, reduction = "umap") + ggtitle(g),
    error = function(e) {
      suppressPackageStartupMessages(library(ks))
      emb <- Embeddings(obj, "umap")[, 1:2]
      w   <- as.numeric(GetAssayData(obj, assay = "RNA", layer = "data")[g, ])
      if (sum(w) == 0) w <- w + 1e-9
      wn  <- w / sum(w) * length(w)
      H   <- ks::Hpi(emb, binned = TRUE)
      dens <- ks::kde(emb, w = wn, H = H, eval.points = emb)$estimate
      df <- data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2], Density = dens)
      df <- df[order(df$Density), ]
      ggplot(df, aes(UMAP_1, UMAP_2, color = Density)) +
        geom_point(size = 0.3) +
        scale_color_viridis_c(name = "Density") +
        ggtitle(g) + theme_classic() +
        theme(plot.title = element_text(hjust = 0.5))
    }
  )
}

plots <- lapply(genes, function(g) density_one(scRNA1, g))
names(plots) <- genes

for (g in genes) {
  ggsave(file.path(RESULTS, paste0("Fig3F_", g, "_density_UMAP.pdf")), plots[[g]], width = 5, height = 4)
  ggsave(file.path(RESULTS, paste0("Fig3F_", g, "_density_UMAP.png")), plots[[g]], width = 5, height = 4, dpi = 300, bg = "white")
}

p_row <- wrap_plots(plots, nrow = 1)
ggsave(file.path(RESULTS, "Fig3F_ABC_transporter_density_UMAP_row.pdf"), p_row, width = 4 * length(genes), height = 4)
ggsave(file.path(RESULTS, "Fig3F_ABC_transporter_density_UMAP_row.png"), p_row, width = 4 * length(genes), height = 4, dpi = 300, bg = "white")
