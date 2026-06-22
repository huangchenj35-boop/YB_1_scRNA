## ============================================================
## Code Ocean stable runner for YB_1_scRNA
## ============================================================
## This runner is designed for a reviewer-facing Code Ocean capsule.
## It avoids hard failure when optional heavy dependencies or external
## intermediate objects are not available.
##
## Default mode:
##   - check and link common /data inputs
##   - run scripts with available dependencies
##   - skip heavy optional steps unless RUN_HEAVY_STEPS=true
##   - write a run log to output/codeocean_run_log.csv
##
## Full mode:
##   RUN_HEAVY_STEPS=true STRICT_CORE=true Rscript run_codeocean.R
## ============================================================

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
set.seed(1234)

as_bool <- function(x, default = FALSE) {
  if (length(x) == 0 || is.na(x) || x == "") {
    return(default)
  }
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

run_heavy_steps <- as_bool(Sys.getenv("RUN_HEAVY_STEPS"), default = FALSE)
strict_core <- as_bool(Sys.getenv("STRICT_CORE"), default = FALSE)
copy_to_results <- as_bool(Sys.getenv("COPY_RESULTS_TO_RESULTS_DIR"), default = TRUE)
allow_no_input <- as_bool(Sys.getenv("ALLOW_NO_INPUT"), default = TRUE)

message("Code Ocean runner settings:")
message("  RUN_HEAVY_STEPS = ", run_heavy_steps)
message("  STRICT_CORE = ", strict_core)
message("  COPY_RESULTS_TO_RESULTS_DIR = ", copy_to_results)
message("  ALLOW_NO_INPUT = ", allow_no_input)

output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

safe_link_or_copy <- function(from, to) {
  if (!file.exists(from) && !dir.exists(from)) {
    return(FALSE)
  }
  if (file.exists(to) || dir.exists(to)) {
    return(TRUE)
  }

  dir.create(dirname(to), showWarnings = FALSE, recursive = TRUE)

  linked <- tryCatch(
    file.symlink(from = from, to = to),
    warning = function(w) FALSE,
    error = function(e) FALSE
  )

  if (isTRUE(linked)) {
    return(TRUE)
  }

  copied <- tryCatch(
    file.copy(from = from, to = to, recursive = TRUE),
    warning = function(w) FALSE,
    error = function(e) FALSE
  )

  isTRUE(copied)
}

prepare_codeocean_inputs <- function() {
  message("\nPreparing Code Ocean input paths...")

  data_roots <- c("/data", "../data", "data")

  for (root in data_roots) {
    gse_path <- file.path(root, "GSE138709")
    if (dir.exists(gse_path) && !dir.exists("GSE138709")) {
      ok <- safe_link_or_copy(gse_path, "GSE138709")
      message("  Linked/copied GSE138709 from ", gse_path, ": ", ok)
      break
    }
  }

  rds_files <- c(
    "scRNA1_preprocessed.rds",
    "scRNA1_annotated.rds",
    "scRNA1_with_SCENIC_AUC.rds",
    "scRNA1_with_DrugPredictions.rds",
    "SCENIC_AUC_matrix.rds",
    "copykat_result.rds",
    "infercnv_result.rds"
  )

  for (rds in rds_files) {
    target <- file.path(output_dir, rds)
    if (file.exists(target)) {
      next
    }

    candidates <- c(
      file.path("/data", rds),
      file.path("/data", "output", rds),
      file.path("../data", rds),
      file.path("../data", "output", rds),
      file.path("data", rds),
      file.path("data", "output", rds)
    )

    candidates <- candidates[file.exists(candidates)]
    if (length(candidates) > 0) {
      ok <- safe_link_or_copy(candidates[[1]], target)
      message("  Linked/copied ", rds, " from ", candidates[[1]], ": ", ok)
    }
  }
}

prepare_codeocean_inputs()

input_found <- any(c(
  dir.exists("GSE138709"),
  dir.exists("Rawcount/filtered_feature_bc_matrix"),
  dir.exists("filtered_feature_bc_matrix"),
  file.exists(file.path(output_dir, "scRNA1_preprocessed.rds")),
  file.exists(file.path(output_dir, "scRNA1_annotated.rds"))
))

if (!input_found) {
  msg <- c(
    "No input data were found.",
    "Place raw GSE138709 10x folders under /data/GSE138709 or provide processed RDS files such as:",
    "  /data/scRNA1_preprocessed.rds",
    "  /data/output/scRNA1_preprocessed.rds",
    "  /data/scRNA1_annotated.rds",
    "The runner exits without error in reviewer-safe mode so that the capsule documents the missing input clearly."
  )
  writeLines(msg, con = file.path(output_dir, "CODEOCEAN_INPUT_REQUIRED.txt"))
  message(paste(msg, collapse = "\n"))
  if (!allow_no_input) {
    stop("No input data were found and ALLOW_NO_INPUT=false.")
  }
  quit(save = "no", status = 0)
}

steps <- list(
  list(script = "00_data_preprocessing_for_FigS1A.R", group = "core", heavy = FALSE,
       packages = c("Seurat", "dplyr", "ggplot2", "patchwork")),
  list(script = "01_Fig1A_FigS1B_cell_annotation.R", group = "core", heavy = FALSE,
       packages = c("Seurat", "dplyr", "ggplot2")),
  list(script = "02_Fig1B_Fig1D_FigS1C_sample_origin_YBX1_composition.R", group = "core", heavy = FALSE,
       packages = c("Seurat", "dplyr", "tidyr", "ggplot2", "scales", "colorspace", "aplot")),
  list(script = "03_Fig1C_Fig1F_Fig3Btop_Fig3Ctop_YBX1_feature_UMAP.R", group = "core", heavy = FALSE,
       packages = c("Seurat", "dplyr", "ggplot2")),
  list(script = "04_Fig1E_FigS2_inferCNV_CopyKAT.R", group = "heavy", heavy = TRUE,
       packages = c("Seurat", "infercnv", "copykat")),
  list(script = "05_Fig1G_Fig3Bbottom_Fig3Cbottom_Fig3D_Fig3E_boxdensity.R", group = "core_or_intermediate", heavy = FALSE,
       packages = c("Seurat", "dplyr", "ggplot2", "ggpubr")),
  list(script = "06_Fig1H_Fig1I_Fig1J_Fig1K_Monocle3_trajectory.R", group = "heavy", heavy = TRUE,
       packages = c("Seurat", "monocle3")),
  list(script = "07_data_drug_sensitivity_prediction_for_Fig3D_Fig3E_FigS3C.R", group = "optional", heavy = FALSE,
       packages = c("Seurat", "dplyr")),
  list(script = "08_data_SCENIC_regulon_inference_for_Fig3_FigS3.R", group = "heavy", heavy = TRUE,
       packages = c("SCENIC", "AUCell", "RcisTarget", "GENIE3")),
  list(script = "09_data_SCENIC_AUC_integration_for_Fig3_FigS3.R", group = "heavy", heavy = TRUE,
       packages = c("Seurat", "AUCell")),
  list(script = "10_Fig3A_SCENIC_CopyKAT_regulon_heatmap.R", group = "heavy", heavy = TRUE,
       packages = c("ComplexHeatmap", "dplyr")),
  list(script = "11_Fig3Ctop_GSVA_YBX1_targets_ssGSEA.R", group = "optional", heavy = FALSE,
       packages = c("Seurat", "GSVA")),
  list(script = "12_Fig3F_FigS3E_to_FigS3H_ABC_transporter.R", group = "core_or_intermediate", heavy = FALSE,
       packages = c("Seurat", "dplyr", "ggplot2")),
  list(script = "13_FigS1D_FigS1E_FigS1F_cholangiocyte_subclustering.R", group = "core_or_intermediate", heavy = FALSE,
       packages = c("Seurat", "dplyr", "ggplot2")),
  list(script = "14_FigS3A_SCENIC_regulon_volcano.R", group = "heavy", heavy = TRUE,
       packages = c("dplyr", "ggplot2")),
  list(script = "15_FigS3B_SCENIC_cisplatin_regulon_heatmap.R", group = "heavy", heavy = TRUE,
       packages = c("ComplexHeatmap")),
  list(script = "16_FigS3C_cisplatin_IC50_UMAP_density.R", group = "optional", heavy = FALSE,
       packages = c("Seurat", "ggplot2")),
  list(script = "17_FigS3D_YBX1_regulon_density_UMAP.R", group = "heavy", heavy = TRUE,
       packages = c("Seurat", "Nebulosa"))
)

run_log <- list()

append_log <- function(script, group, status, reason, start_time, end_time) {
  run_log[[length(run_log) + 1]] <<- data.frame(
    script = script,
    group = group,
    status = status,
    reason = reason,
    start_time = as.character(start_time),
    end_time = as.character(end_time),
    stringsAsFactors = FALSE
  )
}

run_step <- function(step) {
  script <- step$script
  start_time <- Sys.time()

  if (!file.exists(script)) {
    reason <- "script file not found"
    message("\nSKIP: ", script, " -- ", reason)
    append_log(script, step$group, "skipped", reason, start_time, Sys.time())
    return(invisible(FALSE))
  }

  if (isTRUE(step$heavy) && !run_heavy_steps) {
    reason <- "heavy optional step disabled; set RUN_HEAVY_STEPS=true to run"
    message("\nSKIP: ", script, " -- ", reason)
    append_log(script, step$group, "skipped", reason, start_time, Sys.time())
    return(invisible(FALSE))
  }

  missing_pkgs <- step$packages[!vapply(step$packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing_pkgs) > 0) {
    reason <- paste("missing packages:", paste(missing_pkgs, collapse = ", "))
    message("\nSKIP: ", script, " -- ", reason)
    append_log(script, step$group, "skipped", reason, start_time, Sys.time())
    return(invisible(FALSE))
  }

  message("\n============================================================")
  message("Running: ", script)
  message("============================================================")

  ok <- tryCatch(
    {
      source(script, local = .GlobalEnv)
      TRUE
    },
    error = function(e) {
      reason <- conditionMessage(e)
      message("ERROR in ", script, ": ", reason)
      append_log(script, step$group, "failed", reason, start_time, Sys.time())
      FALSE
    }
  )

  if (isTRUE(ok)) {
    append_log(script, step$group, "completed", "", start_time, Sys.time())
  }

  invisible(ok)
}

for (step in steps) {
  run_step(step)
}

run_log_df <- if (length(run_log) > 0) {
  do.call(rbind, run_log)
} else {
  data.frame(
    script = character(),
    group = character(),
    status = character(),
    reason = character(),
    start_time = character(),
    end_time = character(),
    stringsAsFactors = FALSE
  )
}

utils::write.csv(
  run_log_df,
  file = file.path(output_dir, "codeocean_run_log.csv"),
  row.names = FALSE
)

if (copy_to_results && dir.exists("/results")) {
  result_target <- file.path("/results", "output")
  dir.create(result_target, showWarnings = FALSE, recursive = TRUE)
  output_files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
  for (f in output_files) {
    rel <- sub(paste0("^", output_dir, "/?"), "", f)
    target <- file.path(result_target, rel)
    dir.create(dirname(target), showWarnings = FALSE, recursive = TRUE)
    try(file.copy(f, target, overwrite = TRUE), silent = TRUE)
  }
}

message("\nCode Ocean run finished.")
message("Run log: ", file.path(output_dir, "codeocean_run_log.csv"))

failed_core <- run_log_df$status == "failed" & run_log_df$group == "core"
if (strict_core && any(failed_core)) {
  stop("One or more core scripts failed. See output/codeocean_run_log.csv.")
}
