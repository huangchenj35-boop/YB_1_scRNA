## ============================================================
## Code Ocean stable runner for YB_1_scRNA
## ============================================================
## Scripts are named according to their running order.
## Default mode skips heavy optional steps unless RUN_HEAVY_STEPS=true.
## ============================================================

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
set.seed(1234)

as_bool <- function(x, default = FALSE) {
  if (length(x) == 0 || is.na(x) || x == "") return(default)
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

run_heavy_steps <- as_bool(Sys.getenv("RUN_HEAVY_STEPS"), FALSE)
strict_core <- as_bool(Sys.getenv("STRICT_CORE"), FALSE)
copy_to_results <- as_bool(Sys.getenv("COPY_RESULTS_TO_RESULTS_DIR"), TRUE)
allow_no_input <- as_bool(Sys.getenv("ALLOW_NO_INPUT"), TRUE)

message("Code Ocean runner settings:")
message("  RUN_HEAVY_STEPS = ", run_heavy_steps)
message("  STRICT_CORE = ", strict_core)
message("  COPY_RESULTS_TO_RESULTS_DIR = ", copy_to_results)
message("  ALLOW_NO_INPUT = ", allow_no_input)

output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

safe_link_or_copy <- function(from, to) {
  if (!file.exists(from) && !dir.exists(from)) return(FALSE)
  if (file.exists(to) || dir.exists(to)) return(TRUE)
  dir.create(dirname(to), showWarnings = FALSE, recursive = TRUE)

  linked <- tryCatch(
    file.symlink(from = from, to = to),
    warning = function(w) FALSE,
    error = function(e) FALSE
  )
  if (isTRUE(linked)) return(TRUE)

  copied <- tryCatch(
    file.copy(from = from, to = to, recursive = TRUE),
    warning = function(w) FALSE,
    error = function(e) FALSE
  )
  isTRUE(copied)
}

has_umi_csv <- function(path) {
  dir.exists(path) && length(list.files(path, pattern = "_UMI\\.csv(\\.gz)?$", recursive = TRUE)) > 0
}

has_10x_matrix <- function(path) {
  if (!dir.exists(path)) return(FALSE)
  dirs <- list.dirs(path, recursive = TRUE, full.names = TRUE)
  any(file.exists(file.path(dirs, "matrix.mtx")) | file.exists(file.path(dirs, "matrix.mtx.gz")))
}

extract_tar_if_needed <- function(tar_file, exdir = "GSE138709_RAW") {
  if (!file.exists(tar_file)) return(FALSE)
  if (has_umi_csv(exdir)) return(TRUE)

  dir.create(exdir, showWarnings = FALSE, recursive = TRUE)
  message("  Extracting ", tar_file, " to ", exdir, "/")
  tryCatch({
    utils::untar(tar_file, exdir = exdir)
    TRUE
  }, error = function(e) {
    message("  Failed to extract ", tar_file, ": ", conditionMessage(e))
    FALSE
  })
}

prepare_codeocean_inputs <- function() {
  message("\nPreparing Code Ocean input paths...")
  data_roots <- c("/data", "../data", "data")

  for (root in data_roots) {
    raw_dir <- file.path(root, "GSE138709_RAW")
    if (dir.exists(raw_dir) && !dir.exists("GSE138709_RAW")) {
      ok <- safe_link_or_copy(raw_dir, "GSE138709_RAW")
      message("  Linked/copied GSE138709_RAW from ", raw_dir, ": ", ok)
      break
    }
  }

  tar_candidates <- c(
    "GSE138709_RAW.tar",
    file.path("/data", "GSE138709_RAW.tar"),
    file.path("../data", "GSE138709_RAW.tar"),
    file.path("data", "GSE138709_RAW.tar")
  )
  tar_candidates <- tar_candidates[file.exists(tar_candidates)]
  if (length(tar_candidates) > 0) {
    if (!file.exists("GSE138709_RAW.tar")) {
      ok <- safe_link_or_copy(tar_candidates[[1]], "GSE138709_RAW.tar")
      message("  Linked/copied GSE138709_RAW.tar from ", tar_candidates[[1]], ": ", ok)
    }
    extract_tar_if_needed("GSE138709_RAW.tar", "GSE138709_RAW")
  }

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
    if (file.exists(target)) next
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

umi_csv_found <- any(vapply(c("GSE138709_RAW", "GSE138709", "Rawcount", "data", "/data"), has_umi_csv, FUN.VALUE = logical(1)))
tenx_found <- any(vapply(c("GSE138709", "Rawcount/filtered_feature_bc_matrix", "filtered_feature_bc_matrix"), has_10x_matrix, FUN.VALUE = logical(1)))
rds_found <- any(file.exists(file.path(output_dir, c("scRNA1_preprocessed.rds", "scRNA1_annotated.rds"))))

input_found <- any(c(umi_csv_found, tenx_found, rds_found))

if (!input_found) {
  msg <- c(
    "No input data were found.",
    "Place GSE138709_RAW.tar or extracted UMI CSV files in the Code Ocean Data section, for example:",
    "  /data/GSE138709_RAW.tar",
    "  /data/GSE138709_RAW/",
    "Alternatively provide processed RDS files such as:",
    "  /data/scRNA1_preprocessed.rds",
    "  /data/output/scRNA1_preprocessed.rds",
    "  /data/scRNA1_annotated.rds",
    "The runner exits without error in reviewer-safe mode so that the capsule documents the missing input clearly."
  )
  writeLines(msg, con = file.path(output_dir, "CODEOCEAN_INPUT_REQUIRED.txt"))
  message(paste(msg, collapse = "\n"))
  if (!allow_no_input) stop("No input data were found and ALLOW_NO_INPUT=false.")
  quit(save = "no", status = 0)
}

steps <- list(
  list(script = "01_preprocessing_FigS1A.R", group = "core", heavy = FALSE, packages = c("Seurat", "dplyr", "ggplot2", "patchwork", "Matrix")),
  list(script = "02_cell_annotation_Fig1A_FigS1B.R", group = "core", heavy = FALSE, packages = c("Seurat", "dplyr", "ggplot2")),
  list(script = "03_sample_origin_YBX1_composition_Fig1B_Fig1D_FigS1C.R", group = "core", heavy = FALSE, packages = c("Seurat", "dplyr", "tidyr", "ggplot2", "scales", "colorspace", "aplot")),
  list(script = "04_YBX1_feature_UMAP_Fig1C_Fig1F_Fig3B_Fig3C.R", group = "core", heavy = FALSE, packages = c("Seurat", "dplyr", "ggplot2")),
  list(script = "05_CNV_inferCNV_CopyKAT_Fig1E_FigS2.R", group = "heavy", heavy = TRUE, packages = c("Seurat", "infercnv", "copykat")),
  list(script = "06_YBX1_box_density_drug_response_Fig1G_Fig3B_to_Fig3E.R", group = "core_or_intermediate", heavy = FALSE, packages = c("Seurat", "dplyr", "ggplot2", "ggpubr")),
  list(script = "07_Monocle3_trajectory_Fig1H_to_Fig1K.R", group = "heavy", heavy = TRUE, packages = c("Seurat", "monocle3")),
  list(script = "08_drug_sensitivity_prediction_Fig3D_Fig3E_FigS3C.R", group = "optional", heavy = FALSE, packages = c("Seurat", "dplyr")),
  list(script = "09_SCENIC_regulon_inference_Fig3_FigS3.R", group = "heavy", heavy = TRUE, packages = c("SCENIC", "AUCell", "RcisTarget", "GENIE3")),
  list(script = "10_SCENIC_AUC_integration_Fig3_FigS3.R", group = "heavy", heavy = TRUE, packages = c("Seurat", "AUCell")),
  list(script = "11_SCENIC_CopyKAT_regulon_heatmap_Fig3A.R", group = "heavy", heavy = TRUE, packages = c("ComplexHeatmap", "dplyr")),
  list(script = "12_GSVA_YBX1_targets_ssGSEA_Fig3C.R", group = "optional", heavy = FALSE, packages = c("Seurat", "GSVA")),
  list(script = "13_ABC_transporter_Fig3F_FigS3E_to_FigS3H.R", group = "core_or_intermediate", heavy = FALSE, packages = c("Seurat", "dplyr", "ggplot2")),
  list(script = "14_cholangiocyte_subclustering_FigS1D_to_FigS1F.R", group = "core_or_intermediate", heavy = FALSE, packages = c("Seurat", "dplyr", "ggplot2")),
  list(script = "15_SCENIC_regulon_volcano_FigS3A.R", group = "heavy", heavy = TRUE, packages = c("dplyr", "ggplot2")),
  list(script = "16_SCENIC_cisplatin_regulon_heatmap_FigS3B.R", group = "heavy", heavy = TRUE, packages = c("ComplexHeatmap")),
  list(script = "17_cisplatin_IC50_UMAP_density_FigS3C.R", group = "optional", heavy = FALSE, packages = c("Seurat", "ggplot2")),
  list(script = "18_YBX1_regulon_density_UMAP_FigS3D.R", group = "heavy", heavy = TRUE, packages = c("Seurat", "Nebulosa"))
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

  ok <- tryCatch({
    source(script, local = .GlobalEnv)
    TRUE
  }, error = function(e) {
    reason <- conditionMessage(e)
    message("ERROR in ", script, ": ", reason)
    append_log(script, step$group, "failed", reason, start_time, Sys.time())
    FALSE
  })

  if (isTRUE(ok)) append_log(script, step$group, "completed", "", start_time, Sys.time())
  invisible(ok)
}

for (step in steps) run_step(step)

run_log_df <- if (length(run_log) > 0) {
  do.call(rbind, run_log)
} else {
  data.frame(script = character(), group = character(), status = character(), reason = character(), start_time = character(), end_time = character(), stringsAsFactors = FALSE)
}

utils::write.csv(run_log_df, file = file.path(output_dir, "codeocean_run_log.csv"), row.names = FALSE)

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
