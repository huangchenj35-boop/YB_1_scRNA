# YB_1_scRNA

This repository contains the R scripts used for the single-cell RNA-seq analysis of YB-1/YBX1 in intrahepatic cholangiocarcinoma.

The code mainly reproduces the single-cell results shown in Fig. 1, Fig. 3, Fig. S1, Fig. S2, and Fig. S3 of the manuscript. Clinical analyses, immunohistochemistry, and wet-lab validation experiments are not included in this repository.

## Dataset

The analysis was based on the public single-cell RNA-seq dataset **GSE138709**, which includes intrahepatic cholangiocarcinoma tissues and adjacent liver tissues.

Raw count matrices are not included in this repository because of file-size limitations. To run the workflow from raw data, please download the 10x Genomics count matrices from GEO and place the `filtered_feature_bc_matrix` folders under one of the following directories:

```text
GSE138709/
Rawcount/filtered_feature_bc_matrix/
filtered_feature_bc_matrix/
```

Alternatively, users may provide a processed Seurat object at:

```text
output/scRNA1_preprocessed.rds
```

If this file is available, downstream scripts will use it directly.

## Main workflow

The scripts are numbered according to the recommended running order.

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

## Script description

### 1. Preprocessing and integration

```text
00_data_preprocessing_for_FigS1A.R
```

This script reads the single-cell count matrices, creates Seurat objects, applies quality control, integrates samples, performs dimensionality reduction, and generates the cluster UMAP.

The main QC filters are:

```text
nFeature_RNA > 200
nFeature_RNA < 8000
percent.mt < 20
percent.HB < 5
nCount_RNA < 100000
```

Main outputs:

```text
output/scRNA1_preprocessed.rds
output/FigS1A_cluster_UMAP.pdf
output/FigS1A_cluster_UMAP.png
```

### 2. Cell annotation, sample origin, and YBX1 expression

```text
01_Fig1A_FigS1B_cell_annotation.R
02_Fig1B_Fig1D_FigS1C_sample_origin_YBX1_composition.R
03_Fig1C_Fig1F_Fig3Btop_Fig3Ctop_YBX1_feature_UMAP.R
```

These scripts generate the annotated single-cell map, marker-gene plots, tumor/adjacent distribution, cell-type composition, and YBX1 expression plots.

Related panels:

```text
Fig. 1A
Fig. 1B
Fig. 1C
Fig. 1D
Fig. 1F
Fig. S1B
Fig. S1C
Fig. 3B, top
Fig. 3C, top
```

### 3. Malignant-cell inference

```text
04_Fig1E_FigS2_inferCNV_CopyKAT.R
```

This script uses copy-number variation signals to distinguish malignant and non-malignant epithelial cells. It includes inferCNV and CopyKAT-related analyses.

Related panels:

```text
Fig. 1E
Fig. S2
```

### 4. YBX1 level, regulon activity, and predicted drug response

```text
05_Fig1G_Fig3Bbottom_Fig3Cbottom_Fig3D_Fig3E_boxdensity.R
07_data_drug_sensitivity_prediction_for_Fig3D_Fig3E_FigS3C.R
16_FigS3C_cisplatin_IC50_UMAP_density.R
17_FigS3D_YBX1_regulon_density_UMAP.R
```

These scripts compare YBX1 expression, YBX1 regulon activity, YBX1 target-gene scores, and predicted chemotherapy sensitivity.

Related panels:

```text
Fig. 1G
Fig. 3B, bottom
Fig. 3C, bottom
Fig. 3D
Fig. 3E
Fig. S3C
Fig. S3D
```

### 5. Trajectory analysis

```text
06_Fig1H_Fig1I_Fig1J_Fig1K_Monocle3_trajectory.R
```

This script performs Monocle3 trajectory analysis and examines the distribution of malignant-cell status and YBX1 expression along pseudotime.

Related panels:

```text
Fig. 1H
Fig. 1I
Fig. 1J
Fig. 1K
```

### 6. SCENIC regulon analysis

```text
08_data_SCENIC_regulon_inference_for_Fig3_FigS3.R
09_data_SCENIC_AUC_integration_for_Fig3_FigS3.R
10_Fig3A_SCENIC_CopyKAT_regulon_heatmap.R
14_FigS3A_SCENIC_regulon_volcano.R
15_FigS3B_SCENIC_cisplatin_regulon_heatmap.R
```

These scripts perform SCENIC regulon analysis, integrate regulon AUC scores into the Seurat object, and generate regulon heatmaps and volcano plots.

Related panels:

```text
Fig. 3A
Fig. S3A
Fig. S3B
```

### 7. YBX1 target-gene score

```text
11_Fig3Ctop_GSVA_YBX1_targets_ssGSEA.R
```

This script calculates the YBX1 target-gene score using GSVA/ssGSEA and projects the score onto the single-cell embedding.

Related panel:

```text
Fig. 3C, top
```

### 8. ABC transporter genes

```text
12_Fig3F_FigS3E_to_FigS3H_ABC_transporter.R
```

This script examines the expression patterns of selected ABC transporter-related genes.

Genes included:

```text
ABCB1
ABCC1
ABCC2
MVP
```

Related panels:

```text
Fig. 3F
Fig. S3E
Fig. S3F
Fig. S3G
Fig. S3H
```

### 9. Cholangiocyte subclustering

```text
13_FigS1D_FigS1E_FigS1F_cholangiocyte_subclustering.R
```

This script subclusters cholangiocyte-related cells and annotates LPC, cholangiocyte, and malignant-cell populations.

Related panels:

```text
Fig. S1D
Fig. S1E
Fig. S1F
```

## Required R packages

The main R packages used in this workflow include:

```text
Seurat
dplyr
ggplot2
patchwork
Matrix
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
tidyr
readr
stringr
```

SCENIC-related scripts require the corresponding motif ranking databases. These database files should be prepared according to the SCENIC documentation before running the regulon inference step.

## Running the analysis

Run the scripts in numerical order. For example:

```r
source("00_data_preprocessing_for_FigS1A.R")
source("01_Fig1A_FigS1B_cell_annotation.R")
source("02_Fig1B_Fig1D_FigS1C_sample_origin_YBX1_composition.R")
source("03_Fig1C_Fig1F_Fig3Btop_Fig3Ctop_YBX1_feature_UMAP.R")
source("04_Fig1E_FigS2_inferCNV_CopyKAT.R")
source("05_Fig1G_Fig3Bbottom_Fig3Cbottom_Fig3D_Fig3E_boxdensity.R")
source("06_Fig1H_Fig1I_Fig1J_Fig1K_Monocle3_trajectory.R")
source("07_data_drug_sensitivity_prediction_for_Fig3D_Fig3E_FigS3C.R")
source("08_data_SCENIC_regulon_inference_for_Fig3_FigS3.R")
source("09_data_SCENIC_AUC_integration_for_Fig3_FigS3.R")
source("10_Fig3A_SCENIC_CopyKAT_regulon_heatmap.R")
source("11_Fig3Ctop_GSVA_YBX1_targets_ssGSEA.R")
source("12_Fig3F_FigS3E_to_FigS3H_ABC_transporter.R")
source("13_FigS1D_FigS1E_FigS1F_cholangiocyte_subclustering.R")
source("14_FigS3A_SCENIC_regulon_volcano.R")
source("15_FigS3B_SCENIC_cisplatin_regulon_heatmap.R")
source("16_FigS3C_cisplatin_IC50_UMAP_density.R")
source("17_FigS3D_YBX1_regulon_density_UMAP.R")
```

For Code Ocean, these commands can be placed in a `run_all.R` file and executed as:

```bash
Rscript run_all.R
```

## Output files

The scripts write intermediate objects and figures to the `output/` directory.

Common output formats include:

```text
.rds
.pdf
.png
```

Several downstream scripts require intermediate files generated by earlier scripts. If an error reports that an RDS object or regulon AUC matrix is missing, please run the corresponding upstream script first.

## Notes

This repository is meant to document and reproduce the single-cell computational analyses. It does not include raw sequencing matrices, large intermediate RDS files, clinical datasets, IHC images, or wet-lab experimental data.
