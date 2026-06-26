## ============================================================
## Fig S2 — inferCNV genome-wide CNV heatmap (re-plot from final object)
## ------------------------------------------------------------
## Source : 04_FigS2_inferCNV_heatmap.R (infercnv::plot_cnv, ported verbatim).
##          Adapted: object path -> data/ (infercnv_obj.rds kept in the slim
##          data set), output -> results/.
## Needs  : Bioconductor package `infercnv`.
## Output : results/FigS2_inferCNV_heatmap.pdf
## ============================================================

suppressPackageStartupMessages({ library(infercnv) })

if (!exists("RESULTS")) { RESULTS <- Sys.getenv("YB1_RESULTS", unset="../results"); dir.create(RESULTS, recursive=TRUE, showWarnings=FALSE) }
DATA <- Sys.getenv("YB1_DATA", unset = "../data")

## final inferCNV object (run.final.infercnv_obj == infercnv_obj.rds content)
cand <- c(file.path(DATA, "run.final.infercnv_obj"), file.path(DATA, "infercnv_obj.rds"))
infercnv_obj_path <- cand[file.exists(cand)][1]
if (is.na(infercnv_obj_path)) stop("inferCNV object not found in data/ (infercnv_obj.rds).")
infercnv_obj <- readRDS(infercnv_obj_path)

infercnv::plot_cnv(
  infercnv_obj,
  out_dir           = RESULTS,
  output_filename   = "FigS2_inferCNV_heatmap",
  plot_chr_scale    = FALSE,
  cluster_by_groups = TRUE,
  color_safe_pal    = FALSE,
  x.center          = 1,
  x.range           = "auto",
  title             = "inferCNV",
  output_format     = "pdf"
)
message("Fig S2 done: inferCNV heatmap -> ", RESULTS)
