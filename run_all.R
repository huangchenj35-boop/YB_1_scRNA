## ============================================================
## Run all single-cell RNA-seq analysis scripts
## Project: YB_1_scRNA
## ============================================================
## Notes:
##   - Run this file from the repository root.
##   - Raw 10x folders or output/scRNA1_preprocessed.rds should be present
##     before running the full workflow.
##   - SCENIC and drug-response steps require their corresponding input files
##     and databases, as described in the README and script headers.
## ============================================================

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
set.seed(1234)

dir.create("output", showWarnings = FALSE, recursive = TRUE)

scripts <- c(
  "00_data_preprocessing_for_FigS1A.R",
  "01_Fig1A_FigS1B_cell_annotation.R",
  "02_Fig1B_Fig1D_FigS1C_sample_origin_YBX1_composition.R",
  "03_Fig1C_Fig1F_Fig3Btop_Fig3Ctop_YBX1_feature_UMAP.R",
  "04_Fig1E_FigS2_inferCNV_CopyKAT.R",
  "05_Fig1G_Fig3Bbottom_Fig3Cbottom_Fig3D_Fig3E_boxdensity.R",
  "06_Fig1H_Fig1I_Fig1J_Fig1K_Monocle3_trajectory.R",
  "07_data_drug_sensitivity_prediction_for_Fig3D_Fig3E_FigS3C.R",
  "08_data_SCENIC_regulon_inference_for_Fig3_FigS3.R",
  "09_data_SCENIC_AUC_integration_for_Fig3_FigS3.R",
  "10_Fig3A_SCENIC_CopyKAT_regulon_heatmap.R",
  "11_Fig3Ctop_GSVA_YBX1_targets_ssGSEA.R",
  "12_Fig3F_FigS3E_to_FigS3H_ABC_transporter.R",
  "13_FigS1D_FigS1E_FigS1F_cholangiocyte_subclustering.R",
  "14_FigS3A_SCENIC_regulon_volcano.R",
  "15_FigS3B_SCENIC_cisplatin_regulon_heatmap.R",
  "16_FigS3C_cisplatin_IC50_UMAP_density.R",
  "17_FigS3D_YBX1_regulon_density_UMAP.R"
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
