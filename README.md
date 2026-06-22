# YB_1_scRNA

R scripts for the single-cell RNA-seq analysis of YB-1/YBX1 in intrahepatic cholangiocarcinoma.

This repository covers the computational analyses used for the single-cell figures in the manuscript, mainly Fig. 1, Fig. 3, Fig. S1, Fig. S2, and Fig. S3. Clinical analyses, IHC images, wet-lab validation data, and tabulated experimental results are not included.

## Data

The analysis uses the public dataset **GSE138709**, including 5 intrahepatic cholangiocarcinoma tumor samples and 3 adjacent liver samples.

The GEO supplementary data are provided as processed UMI count matrices in `GSE138709_RAW.tar`. The archive contains eight `*_UMI.csv.gz` files:

```text
GSM4116579_ICC_18_Adjacent_UMI.csv.gz
GSM4116580_ICC_18_Tumor_UMI.csv.gz
GSM4116581_ICC_20_Tumor_UMI.csv.gz
GSM4116582_ICC_23_Adjacent_UMI.csv.gz
GSM4116583_ICC_23_Tumor_UMI.csv.gz
GSM4116584_ICC_24_Tumor1_UMI.csv.gz
GSM4116585_ICC_24_Tumor2_UMI.csv.gz
GSM4116586_ICC_25_Adjacent_UMI.csv.gz
```

Raw count matrices and large intermediate RDS files are not stored in this repository. For Code Ocean, upload `GSE138709_RAW.tar` to the Data section as:

```text
/data/GSE138709_RAW.tar
```

Alternatively, upload the extracted UMI matrices as:

```text
/data/GSE138709_RAW/*_UMI.csv.gz
```

The preprocessing script can also read standard 10x Genomics `filtered_feature_bc_matrix` folders if they are provided under one of these paths:

```text
GSE138709/
Rawcount/filtered_feature_bc_matrix/
filtered_feature_bc_matrix/
```

If a processed Seurat object is already available, place it here:

```text
output/scRNA1_preprocessed.rds
```

For Code Ocean, the following Data paths are also recognized by `run_codeocean.R`:

```text
/data/scRNA1_preprocessed.rds
/data/output/scRNA1_preprocessed.rds
/data/scRNA1_annotated.rds
/data/output/scRNA1_annotated.rds
```

Downstream scripts use these objects when present.

## Script order

Run the scripts in numerical order:

```text
00_data_preprocessing_for_FigS1A.R
01_Fig1A_FigS1B_cell_annotation.R
02_Fig1B_Fig1D_FigS1C_sample_origin_YBX1_composition.R
03_Fig1C_Fig1F_Fig3Btop_Fig3Ctop_YBX1_feature_UMAP.R
04_Fig1E_FigS2_inferCNV_CopyKAT.R
05_Fig1G_Fig3Bbottom_Fig3Cbottom_Fig3D_Fig3E_boxdensity.R
06_Fig1H_Fig1I_Fig1J_Fig1K_Monocle3_trajectory.R
07_data_drug_sensitivity_prediction_for_Fig3D_Fig3E_FigS3C.R
08_data_SCENIC_regulon_inference_for_Fig3_FigS3.R
09_data_SCENIC_AUC_integration_for_Fig3_FigS3.R
10_Fig3A_SCENIC_CopyKAT_regulon_heatmap.R
11_Fig3Ctop_GSVA_YBX1_targets_ssGSEA.R
12_Fig3F_FigS3E_to_FigS3H_ABC_transporter.R
13_FigS1D_FigS1E_FigS1F_cholangiocyte_subclustering.R
14_FigS3A_SCENIC_regulon_volcano.R
15_FigS3B_SCENIC_cisplatin_regulon_heatmap.R
16_FigS3C_cisplatin_IC50_UMAP_density.R
17_FigS3D_YBX1_regulon_density_UMAP.R
```

Some scripts generate figures directly. Others prepare intermediate Seurat, SCENIC, or drug-response objects used by later scripts.

## Script map

| Script | Main purpose | Related panels |
|---|---|---|
| `00_data_preprocessing_for_FigS1A.R` | Read GEO UMI CSV matrices or 10x folders; Seurat preprocessing, QC, integration, clustering, UMAP/t-SNE | Fig. S1A |
| `01_Fig1A_FigS1B_cell_annotation.R` | Cell-type annotation and marker validation | Fig. 1A; Fig. S1B |
| `02_Fig1B_Fig1D_FigS1C_sample_origin_YBX1_composition.R` | Tumor/adjacent distribution, YBX1 expression comparison, cell-type composition | Fig. 1B; Fig. 1D; Fig. S1C |
| `03_Fig1C_Fig1F_Fig3Btop_Fig3Ctop_YBX1_feature_UMAP.R` | YBX1 expression and YBX1-related score projection | Fig. 1C; Fig. 1F; Fig. 3B top; Fig. 3C top |
| `04_Fig1E_FigS2_inferCNV_CopyKAT.R` | CNV-based malignant-cell inference | Fig. 1E; Fig. S2 |
| `05_Fig1G_Fig3Bbottom_Fig3Cbottom_Fig3D_Fig3E_boxdensity.R` | YBX1 expression, regulon activity, target score, and drug-response comparison | Fig. 1G; Fig. 3B bottom; Fig. 3C bottom; Fig. 3D; Fig. 3E |
| `06_Fig1H_Fig1I_Fig1J_Fig1K_Monocle3_trajectory.R` | Monocle3 trajectory analysis | Fig. 1H-K |
| `07_data_drug_sensitivity_prediction_for_Fig3D_Fig3E_FigS3C.R` | Predicted cisplatin/gemcitabine sensitivity data | Fig. 3D; Fig. 3E; Fig. S3C |
| `08_data_SCENIC_regulon_inference_for_Fig3_FigS3.R` | SCENIC regulon inference | Fig. 3; Fig. S3 |
| `09_data_SCENIC_AUC_integration_for_Fig3_FigS3.R` | Integration of regulon AUC scores into Seurat metadata | Fig. 3; Fig. S3 |
| `10_Fig3A_SCENIC_CopyKAT_regulon_heatmap.R` | Regulon heatmap by malignant-cell status | Fig. 3A |
| `11_Fig3Ctop_GSVA_YBX1_targets_ssGSEA.R` | GSVA/ssGSEA YBX1 target-gene score | Fig. 3C top |
| `12_Fig3F_FigS3E_to_FigS3H_ABC_transporter.R` | ABC transporter expression analysis | Fig. 3F; Fig. S3E-H |
| `13_FigS1D_FigS1E_FigS1F_cholangiocyte_subclustering.R` | Cholangiocyte-related subclustering | Fig. S1D-F |
| `14_FigS3A_SCENIC_regulon_volcano.R` | Regulon volcano plot | Fig. S3A |
| `15_FigS3B_SCENIC_cisplatin_regulon_heatmap.R` | Cisplatin-related regulon heatmap | Fig. S3B |
| `16_FigS3C_cisplatin_IC50_UMAP_density.R` | Predicted cisplatin lnIC50 density UMAP | Fig. S3C |
| `17_FigS3D_YBX1_regulon_density_UMAP.R` | YBX1 regulon AUC density UMAP | Fig. S3D |

## Preprocessing settings

The preprocessing script applies the following QC filters:

```text
nFeature_RNA > 200
nFeature_RNA < 8000
percent.mt < 20
percent.HB < 5
nCount_RNA < 100000
```

Main preprocessing outputs:

```text
output/scRNA1_preprocessed.rds
output/FigS1A_cluster_UMAP.pdf
output/FigS1A_cluster_UMAP.png
```

## R packages

Main packages used in the scripts:

```text
Seurat
dplyr
tidyr
readr
stringr
ggplot2
patchwork
Matrix
scales
colorspace
aplot
Nebulosa
monocle3
infercnv
copykat
SCENIC
AUCell
RcisTarget
GENIE3
GSVA
ComplexHeatmap
pheatmap
ggpubr
cowplot
```

SCENIC-related scripts also require the corresponding motif-ranking database files. These database files are not included and should be prepared following the SCENIC workflow.

## Running the analysis

For Code Ocean review mode, use the stable runner:

```bash
Rscript run_codeocean.R
```

This runner detects `/data/GSE138709_RAW.tar`, extracts the UMI CSV matrices, runs available core scripts, skips heavy optional steps when dependencies are unavailable, and writes:

```text
output/codeocean_run_log.csv
```

For a strict full workflow in a fully prepared local environment, run:

```bash
Rscript run_all.R
```

The individual scripts can also be run in R with `source()` in numerical order.

## Output

The scripts generate `.rds`, `.pdf`, `.png`, and `.csv` files. Most outputs are written to `output/`; a few plotting scripts write files to the working directory because they follow the original figure-export code.

Downstream scripts often require objects generated by earlier scripts. If a script reports a missing Seurat object, regulon AUC matrix, or predicted drug-response object, run the upstream script first or place the required object in the expected path.

## Repository scope

This repository documents the single-cell computational analysis only. It does not include raw sequencing matrices, large intermediate RDS files, clinical datasets, IHC images, wet-lab experimental data, or manual figure assembly files.
