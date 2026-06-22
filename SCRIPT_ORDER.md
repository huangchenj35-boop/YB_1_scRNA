# Ordered script names

The recommended entry scripts are named according to the running order. These files live in the repository root and are used by `run_all.R` and `run_codeocean.R`.

```text
00_install_core_packages.R
01_preprocessing_FigS1A.R
02_cell_annotation_Fig1A_FigS1B.R
03_sample_origin_YBX1_composition_Fig1B_Fig1D_FigS1C.R
04_YBX1_feature_UMAP_Fig1C_Fig1F_Fig3B_Fig3C.R
05_CNV_inferCNV_CopyKAT_Fig1E_FigS2.R
06_YBX1_box_density_drug_response_Fig1G_Fig3B_to_Fig3E.R
07_Monocle3_trajectory_Fig1H_to_Fig1K.R
08_drug_sensitivity_prediction_Fig3D_Fig3E_FigS3C.R
09_SCENIC_regulon_inference_Fig3_FigS3.R
10_SCENIC_AUC_integration_Fig3_FigS3.R
11_SCENIC_CopyKAT_regulon_heatmap_Fig3A.R
12_GSVA_YBX1_targets_ssGSEA_Fig3C.R
13_ABC_transporter_Fig3F_FigS3E_to_FigS3H.R
14_cholangiocyte_subclustering_FigS1D_to_FigS1F.R
15_SCENIC_regulon_volcano_FigS3A.R
16_SCENIC_cisplatin_regulon_heatmap_FigS3B.R
17_cisplatin_IC50_UMAP_density_FigS3C.R
18_YBX1_regulon_density_UMAP_FigS3D.R
```

The original implementation scripts have been backed up under:

```text
legacy_scripts/
```

Each ordered script calls its corresponding implementation script in `legacy_scripts/` with `source()`.
