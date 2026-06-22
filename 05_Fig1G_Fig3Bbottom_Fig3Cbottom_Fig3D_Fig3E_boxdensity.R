## ============================================================
## YB-1 Metrics Box-Density Comparisons
## Figures: Fig 1G, Fig 3B bottom, Fig 3C bottom, Fig 3D, Fig 3E
## Dataset: GSE138709 (5 tumor + 3 adjacent iCCA samples)
## ============================================================
## Prerequisites:
##   - scRNA1 with group_copykat metadata
##   - Cisplatin_pred_group / Gemcitabine_pred_group from script 07
##   - YBX1 regulon AUC from script 09
##   - GSVA_YBX1_targets_ssgsea from script 11
## Output:
##   - BoxDensity_YBX1_by_group_copykat_AllCells_group_copykat.pdf (Fig 1G)
##   - BoxDensity_<YBX1_AUC>_by_group_copykat_*.pdf                (Fig 3B bottom)
##   - BoxDensity_GSVA_YBX1_targets_ssgsea_by_group_copykat_*.pdf  (Fig 3C bottom)
##   - BoxDensity_*_by_Cisplatin_pred_group_AllCells.pdf           (Fig 3D)
##   - BoxDensity_*_by_Gemcitabine_pred_group_AllCells.pdf         (Fig 3E)
## Notes:
##   - This script only uses existing grouping columns. It does not recalculate
##     or hard-code drug cutoffs.
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
dir.create("output", showWarnings = FALSE, recursive = TRUE)

if (!exists("scRNA1") || !inherits(scRNA1, "Seurat")) {
  candidates <- c(
    "output/scRNA1_with_DrugPredictions.rds",
    "output/scRNAauc.rds",
    "output/scRNA1_annotated.rds",
    "output/scRNA1_preprocessed.rds"
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) {
    stop("scRNA1 was not found in memory and no expected RDS file exists in output/.")
  }
  scRNA1 <- readRDS(hit)
}

safe <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

find_ybx1_gene <- function(obj) {
  assays <- names(obj@assays)
  for (assay in assays) {
    rn <- tryCatch(rownames(obj@assays[[assay]]), error = function(e) character(0))
    if ("YBX1" %in% rn) return(list(assay = assay, feature = "YBX1"))
  }
  stop("YBX1 gene was not found in any assay.")
}

tf_root <- function(x) {
  y <- sub("\\s*\\(\\d+g\\)$", "", x)
  y <- sub("_[0-9]+g$", "", y, ignore.case = TRUE)
  y <- sub("_extended$", "", y, ignore.case = TRUE)
  y
}

load_auc_matrix <- function() {
  if (exists("auc_mtx")) return(auc_mtx)
  if (file.exists("output/auc_mtx_regulon_by_cell.rds")) {
    return(readRDS("output/auc_mtx_regulon_by_cell.rds"))
  }
  NULL
}

ensure_ybx1_auc <- function(obj) {
  mdn <- colnames(obj@meta.data)
  preferred <- intersect(c("YBX1_extended_907g", "YBX1_extended_907G"), mdn)
  if (length(preferred)) return(list(obj = obj, col = preferred[1]))

  md_hits <- grep("^YBX1.*(regulon|AUC|extended|[0-9]+g)", mdn,
                  ignore.case = TRUE, value = TRUE)
  md_hits <- md_hits[tf_root(md_hits) == "YBX1"]
  if (length(md_hits)) return(list(obj = obj, col = md_hits[1]))

  auc <- load_auc_matrix()
  if (!is.null(auc)) {
    common <- intersect(colnames(obj), colnames(auc))
    auc_hits <- rownames(auc)[tf_root(rownames(auc)) == "YBX1"]
    if (length(common) && length(auc_hits)) {
      preferred_auc <- auc_hits[auc_hits == "YBX1_extended_907g"]
      pick <- if (length(preferred_auc)) preferred_auc[1] else auc_hits[1]
      vals <- rep(NA_real_, ncol(obj))
      names(vals) <- colnames(obj)
      vals[common] <- as.numeric(auc[pick, common])
      col <- make.names(pick)
      obj[[col]] <- vals
      return(list(obj = obj, col = col))
    }
  }

  warning("YBX1 regulon AUC was not found; Fig 3B bottom will be skipped.")
  list(obj = obj, col = NA_character_)
}

ensure_gsva <- function(obj) {
  mdn <- colnames(obj@meta.data)
  hits <- grep("^GSVA[_-]?YBX1.*targets.*ssgsea$", mdn,
               ignore.case = TRUE, value = TRUE)
  if (length(hits)) return(list(obj = obj, col = hits[1]))

  warning("GSVA_YBX1_targets_ssgsea was not found; Fig 3C bottom will be skipped.")
  list(obj = obj, col = NA_character_)
}

plot_density_box <- function(obj, feature, feature_label, group_col, group_levels,
                             colors, tag) {
  if (is.na(feature) || !nzchar(feature)) return(invisible(NULL))
  if (!group_col %in% colnames(obj@meta.data)) {
    message("Skipping ", feature_label, " by ", group_col, ": grouping column is missing.")
    return(invisible(NULL))
  }

  if (feature %in% colnames(obj@meta.data)) {
    values <- as.numeric(obj@meta.data[[feature]])
  } else {
    values <- as.numeric(FetchData(obj, vars = feature)[, 1])
  }

  df <- data.frame(
    group = factor(as.character(obj@meta.data[[group_col]]), levels = group_levels),
    value = values,
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$value) & !is.na(df$group), , drop = FALSE]
  if (!nrow(df)) {
    message("Skipping ", feature_label, " by ", group_col, ": no valid values.")
    return(invisible(NULL))
  }

  present_levels <- levels(droplevels(df$group))
  if (length(present_levels) == 2) {
    p_val <- tryCatch(stats::wilcox.test(value ~ group, data = df)$p.value,
                      error = function(e) NA_real_)
  } else if (length(present_levels) > 2) {
    p_val <- tryCatch(stats::kruskal.test(value ~ group, data = df)$p.value,
                      error = function(e) NA_real_)
  } else {
    p_val <- NA_real_
  }
  p_lab <- if (is.finite(p_val)) {
    if (p_val < 0.001) "P < 0.001" else paste0("P = ", signif(p_val, 3))
  } else {
    "P = NA"
  }

  color_use <- colors[group_levels]
  color_use <- color_use[!is.na(color_use)]

  p_top <- ggplot(df, aes(x = value, color = group, fill = group)) +
    geom_density(alpha = 0.25, linewidth = 0.6) +
    geom_rug(alpha = 0.25) +
    scale_color_manual(values = color_use, drop = FALSE) +
    scale_fill_manual(values = alpha(color_use, 0.35), drop = FALSE) +
    labs(x = feature_label, y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    )

  p_bottom <- ggplot(df, aes(x = group, y = value, fill = group, color = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.35) +
    scale_fill_manual(values = color_use, drop = FALSE) +
    scale_color_manual(values = color_use, drop = FALSE) +
    annotate("text", x = 1.5, y = max(df$value, na.rm = TRUE),
             label = p_lab, angle = 270, fontface = "bold", size = 3.4) +
    coord_flip() +
    labs(x = NULL, y = feature_label) +
    theme_void(base_size = 11) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(color = "black")
    )

  p_all <- p_top / p_bottom + patchwork::plot_layout(heights = c(2, 1))
  base <- file.path(
    "output",
    paste0("BoxDensity_", safe(feature_label), "_by_", group_col, "_", tag)
  )
  ggsave(paste0(base, ".pdf"), p_all, width = 6, height = 3.6)
  ggsave(paste0(base, ".png"), p_all, width = 6, height = 3.6, dpi = 300, bg = "white")
  invisible(p_all)
}

ybx1_info <- find_ybx1_gene(scRNA1)
DefaultAssay(scRNA1) <- ybx1_info$assay
ybx1_feature <- ybx1_info$feature

auc_info <- ensure_ybx1_auc(scRNA1)
scRNA1 <- auc_info$obj
ybx1_auc_col <- auc_info$col

gsva_info <- ensure_gsva(scRNA1)
scRNA1 <- gsva_info$obj
gsva_col <- gsva_info$col

features_to_plot <- list(
  list(feature = ybx1_feature, label = "YBX1"),
  list(feature = ybx1_auc_col, label = ybx1_auc_col),
  list(feature = gsva_col, label = gsva_col)
)

copykat_levels <- c("aneuploid", "diploid")
copykat_colors <- c("aneuploid" = "#F29C64", "diploid" = "#5AAFD6")

drug_levels <- c("Predicted sensitive", "Predicted resistant")
drug_colors <- c("Predicted sensitive" = "#377EB8", "Predicted resistant" = "#E41A1C")

if (!"group_copykat" %in% colnames(scRNA1@meta.data)) {
  stop("group_copykat metadata is missing. Run 04_Fig1E_FigS2_inferCNV_CopyKAT.R first.")
}

for (item in features_to_plot) {
  plot_density_box(scRNA1, item$feature, item$label, "group_copykat",
                   copykat_levels, copykat_colors, "AllCells_group_copykat")
}

for (item in features_to_plot) {
  plot_density_box(scRNA1, item$feature, item$label, "Cisplatin_pred_group",
                   drug_levels, drug_colors, "AllCells")
}

for (item in features_to_plot) {
  plot_density_box(scRNA1, item$feature, item$label, "Gemcitabine_pred_group",
                   drug_levels, drug_colors, "AllCells")
}

saveRDS(scRNA1, "output/scRNA1_boxdensity_inputs.rds")
message("Done: Fig 1G, Fig 3B bottom, Fig 3C bottom, Fig 3D, and Fig 3E box-density plots saved.")
