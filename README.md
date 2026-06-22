# YB_1_scRNA
Code for YB-1 scRNA analyses
PDF figure order reviewed:

1. Figure 1, page 45: single-cell landscape and trajectory panels A-K
2. Figure 2, page 46: clinical expression, IHC, and survival panels A-I
3. Figure 3, page 47: YB-1 regulon, chemotherapy prediction, ABC transporter panels A-F
4. Figure 4, page 48 through Figure 7, page 51: experimental validation panels
5. Figure S1, page 52; Figure S2, page 53; Figure S3, page 54
6. Figure S4 through Figure S7, pages 55-58
7. Tables 1-2, pages 59-61

This folder contains the scRNA-seq/single-cell code found in `Code_Final` and
`Code_Final_V2`, reordered to match the article figure order. Figures 2, 4-7,
S4-S7, and Tables 1-2 are clinical, wet-lab, or table outputs; no matching R
scripts were found in the supplied code folders.

For the panel-by-panel audit table, see `FIGURE_PANEL_INDEX.csv`.
For the biological label/grouping checklist, see `CHECKLIST.md`.

## Script order and panel mapping

| V3 script | PDF panel(s) | Content checked against PDF |
|---|---|---|
| `00_data_preprocessing_for_FigS1A.R` | Fig S1A | Integrated Seurat object and cluster UMAP source |
| `01_Fig1A_FigS1B_cell_annotation.R` | Fig 1A, Fig S1B | Article-standard cell-type labels and marker violin plots |
| `02_Fig1B_Fig1D_FigS1C_sample_origin_YBX1_composition.R` | Fig 1B, Fig 1D, Fig S1C | Tumor/adjacent UMAP, YB-1 expression comparison, cell-type proportions |
| `03_Fig1C_Fig1F_Fig3Btop_Fig3Ctop_YBX1_feature_UMAP.R` | Fig 1C, Fig 1F, Fig 3B top, Fig 3C top | YB-1 density, YB-1 expression, YB-1 regulon AUC, and GSVA YB-1 target score UMAPs |
| `04_Fig1E_FigS2_inferCNV_CopyKAT.R` | Fig 1E, Fig S2 | CopyKAT diploid/aneuploid UMAP and inferCNV CNV heatmap |
| `05_Fig1G_Fig3Bbottom_Fig3Cbottom_Fig3D_Fig3E_boxdensity.R` | Fig 1G, Fig 3B bottom, Fig 3C bottom, Fig 3D, Fig 3E | YB-1, regulon AUC, GSVA score, cisplatin, and gemcitabine box-density comparisons |
| `06_Fig1H_Fig1I_Fig1J_Fig1K_Monocle3_trajectory.R` | Fig 1H-K | Pseudotime, CNV status, trajectory state, and YB-1 along pseudotime |
| `07_data_drug_sensitivity_prediction_for_Fig3D_Fig3E_FigS3C.R` | Data source | Cisplatin/gemcitabine predicted sensitivity groups and IC50 metadata |
| `08_data_SCENIC_regulon_inference_for_Fig3_FigS3.R` | Data source | SCENIC regulon inference |
| `09_data_SCENIC_AUC_integration_for_Fig3_FigS3.R` | Data source | Adds SCENIC AUC scores, including YB-1 regulon AUC, to Seurat metadata |
| `10_Fig3A_SCENIC_CopyKAT_regulon_heatmap.R` | Fig 3A | Aneuploid/diploid regulon activity heatmap |
| `11_Fig3Ctop_GSVA_YBX1_targets_ssGSEA.R` | Fig 3C top | GSVA/ssGSEA YB-1 target score calculation and UMAP |
| `12_Fig3F_FigS3E_to_FigS3H_ABC_transporter.R` | Fig 3F, Fig S3E-H | ABCB1/ABCC1/ABCC2/MVP density UMAPs and YB-1/ABC heatmaps |
| `13_FigS1D_FigS1E_FigS1F_cholangiocyte_subclustering.R` | Fig S1D-F | Cholangiocyte subclusters, LPC/cholangiocyte/malignant labels, subtype markers |
| `14_FigS3A_SCENIC_regulon_volcano.R` | Fig S3A | Cisplatin-sensitive vs resistant regulon volcano |
| `15_FigS3B_SCENIC_cisplatin_regulon_heatmap.R` | Fig S3B | Same cisplatin-sensitive vs resistant grouping as Fig S3A; regulon AUC heatmap |
| `16_FigS3C_cisplatin_IC50_UMAP_density.R` | Fig S3C | Predicted cisplatin lnIC50 density UMAP |
| `17_FigS3D_YBX1_regulon_density_UMAP.R` | Fig S3D | YB-1 regulon AUC density UMAP |
