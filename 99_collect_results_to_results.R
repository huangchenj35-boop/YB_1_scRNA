## ============================================================
## Collect Code Ocean results
## ============================================================
## This script copies generated result files into /results for Code Ocean.
## It does not rename files and does not change result formats.
##
## It can be run after any analysis command:
##   Rscript 99_collect_results_to_results.R
## ============================================================

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

if (!dir.exists("/results")) {
  message("/results does not exist. Nothing to copy. This is expected outside Code Ocean.")
  quit(save = "no", status = 0)
}

copy_one <- function(from, to) {
  dir.create(dirname(to), showWarnings = FALSE, recursive = TRUE)
  ok <- tryCatch(
    file.copy(from, to, overwrite = TRUE),
    warning = function(w) FALSE,
    error = function(e) FALSE
  )
  isTRUE(ok)
}

copied <- character(0)

if (dir.exists("output")) {
  output_files <- list.files("output", recursive = TRUE, full.names = TRUE, all.files = FALSE)
  for (f in output_files) {
    if (!file.info(f)$isdir) {
      rel <- sub("^output/?", "", f)
      target <- file.path("/results", "output", rel)
      if (copy_one(f, target)) copied <- c(copied, target)
    }
  }
}

root_patterns <- c("pdf", "png", "csv", "tsv", "txt", "rds")
root_files <- list.files(".", pattern = paste0("\\.(", paste(root_patterns, collapse = "|"), ")$"), full.names = TRUE, ignore.case = TRUE)
root_files <- root_files[!grepl("^\\./output/", root_files)]
root_files <- root_files[!grepl("^\\./legacy_scripts/", root_files)]
root_files <- root_files[!grepl("^\\./GSE138709", root_files)]

for (f in root_files) {
  target <- file.path("/results", basename(f))
  if (copy_one(f, target)) copied <- c(copied, target)
}

manifest <- data.frame(file = copied, stringsAsFactors = FALSE)
utils::write.csv(manifest, file = "/results/result_manifest.csv", row.names = FALSE)

message("Copied ", length(copied), " result files to /results.")
message("Manifest: /results/result_manifest.csv")
