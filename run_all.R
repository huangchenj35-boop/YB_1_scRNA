## ============================================================
## Run all single-cell RNA-seq analysis scripts
## Project: YB_1_scRNA
## ============================================================
## Notes:
##   - Run this file from the repository root.
##   - Raw UMI matrices, 10x folders, or output/scRNA1_preprocessed.rds
##     should be present before running the full workflow.
##   - Execution follows data dependency order. Output filenames and figure
##     formats are defined inside each analysis script and are not changed here.
## ============================================================

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
set.seed(1234)

dir.create("output", showWarnings = FALSE, recursive = TRUE)

scripts <- c(
  "01_preprocessing_FigS1A.R",
  "02_cell_annotation_Fig1A_FigS1B.R",
  "03_sample_origin_YBX1_composition_Fig1B_Fig1D_FigS1C.R",
  "04_YBX1_feature_UMAP_Fig1C_Fig1F_Fig3B_Fig3C.R",
  "05_CNV_inferCNV_CopyKAT_Fig1E_FigS2.R",
  "07_Monocle3_trajectory_Fig1H_to_Fig1K.R",
  "08_drug_sensitivity_prediction_Fig3D_Fig3E_FigS3C.R",
  "09_SCENIC_regulon_inference_Fig3_FigS3.R",
  "10_SCENIC_AUC_integration_Fig3_FigS3.R",
  "12_GSVA_YBX1_targets_ssGSEA_Fig3C.R",
  "06_YBX1_box_density_drug_response_Fig1G_Fig3B_to_Fig3E.R",
  "11_SCENIC_CopyKAT_regulon_heatmap_Fig3A.R",
  "13_ABC_transporter_Fig3F_FigS3E_to_FigS3H.R",
  "14_cholangiocyte_subclustering_FigS1D_to_FigS1F.R",
  "15_SCENIC_regulon_volcano_FigS3A.R",
  "16_SCENIC_cisplatin_regulon_heatmap_FigS3B.R",
  "17_cisplatin_IC50_UMAP_density_FigS3C.R",
  "18_YBX1_regulon_density_UMAP_FigS3D.R"
)

missing_scripts <- scripts[!file.exists(scripts)]
if (length(missing_scripts) > 0) {
  stop(
    "The following scripts are missing from the repository root: ",
    paste(missing_scripts, collapse = ", ")
  )
}

for (script in scripts) {
  message("\n============================================================")
  message("Running: ", script)
  message("============================================================")
  source(script, local = .GlobalEnv)
}

message("\nAll listed scripts have finished.")
