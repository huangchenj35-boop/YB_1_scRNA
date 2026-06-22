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

For Code Ocean, upload `GSE138709_RAW.tar` to the Data section as:

```text
/data/GSE138709_RAW.tar
```

Alternatively, upload the extracted UMI matrices as:

```text
/data/GSE138709_RAW/
```

If a processed Seurat object is already available, the following paths are recognized:

```text
output/scRNA1_preprocessed.rds
/data/scRNA1_preprocessed.rds
/data/output/scRNA1_preprocessed.rds
/data/scRNA1_annotated.rds
/data/output/scRNA1_annotated.rds
```

## Installation and run order

For Code Ocean review mode, use the following order in a clean capsule:

```bash
bash codeocean_system_deps.sh
Rscript 00_install_core_packages.R
Rscript packages.R
Rscript run_codeocean.R
```

If system libraries are already available, skip the first line:

```bash
Rscript 00_install_core_packages.R
Rscript packages.R
Rscript run_codeocean.R
```

The main Code Ocean run command can remain:

```bash
Rscript run_codeocean.R
```

The installation details are documented in `INSTALL.md`.

## Ordered script entry points

Scripts are named according to the running order. These ordered files are the recommended entry points for local use and Code Ocean.

```text
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

The implementation scripts are backed up under `legacy_scripts/`. Each ordered entry script calls the corresponding implementation script with `source()`.

## Script map

| Order | Ordered script | Main purpose | Related panels |
|---|---|---|---|
| 01 | `01_preprocessing_FigS1A.R` | Read GEO UMI CSV matrices or 10x folders; Seurat preprocessing, QC, integration, clustering, UMAP/t-SNE | Fig. S1A |
| 02 | `02_cell_annotation_Fig1A_FigS1B.R` | Cell-type annotation and marker validation | Fig. 1A; Fig. S1B |
| 03 | `03_sample_origin_YBX1_composition_Fig1B_Fig1D_FigS1C.R` | Tumor/adjacent distribution, YBX1 expression comparison, cell-type composition | Fig. 1B; Fig. 1D; Fig. S1C |
| 04 | `04_YBX1_feature_UMAP_Fig1C_Fig1F_Fig3B_Fig3C.R` | YBX1 expression and YBX1-related score projection | Fig. 1C; Fig. 1F; Fig. 3B top; Fig. 3C top |
| 05 | `05_CNV_inferCNV_CopyKAT_Fig1E_FigS2.R` | CNV-based malignant-cell inference | Fig. 1E; Fig. S2 |
| 06 | `06_YBX1_box_density_drug_response_Fig1G_Fig3B_to_Fig3E.R` | YBX1 expression, regulon activity, target score, and drug-response comparison | Fig. 1G; Fig. 3B bottom; Fig. 3C bottom; Fig. 3D; Fig. 3E |
| 07 | `07_Monocle3_trajectory_Fig1H_to_Fig1K.R` | Monocle3 trajectory analysis | Fig. 1H-K |
| 08 | `08_drug_sensitivity_prediction_Fig3D_Fig3E_FigS3C.R` | Predicted cisplatin/gemcitabine sensitivity data | Fig. 3D; Fig. 3E; Fig. S3C |
| 09 | `09_SCENIC_regulon_inference_Fig3_FigS3.R` | SCENIC regulon inference | Fig. 3; Fig. S3 |
| 10 | `10_SCENIC_AUC_integration_Fig3_FigS3.R` | Integration of regulon AUC scores into Seurat metadata | Fig. 3; Fig. S3 |
| 11 | `11_SCENIC_CopyKAT_regulon_heatmap_Fig3A.R` | Regulon heatmap by malignant-cell status | Fig. 3A |
| 12 | `12_GSVA_YBX1_targets_ssGSEA_Fig3C.R` | GSVA/ssGSEA YBX1 target-gene score | Fig. 3C top |
| 13 | `13_ABC_transporter_Fig3F_FigS3E_to_FigS3H.R` | ABC transporter expression analysis | Fig. 3F; Fig. S3E-H |
| 14 | `14_cholangiocyte_subclustering_FigS1D_to_FigS1F.R` | Cholangiocyte-related subclustering | Fig. S1D-F |
| 15 | `15_SCENIC_regulon_volcano_FigS3A.R` | Regulon volcano plot | Fig. S3A |
| 16 | `16_SCENIC_cisplatin_regulon_heatmap_FigS3B.R` | Cisplatin-related regulon heatmap | Fig. S3B |
| 17 | `17_cisplatin_IC50_UMAP_density_FigS3C.R` | Predicted cisplatin lnIC50 density UMAP | Fig. S3C |
| 18 | `18_YBX1_regulon_density_UMAP_FigS3D.R` | YBX1 regulon AUC density UMAP | Fig. S3D |

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

Core packages installed by `00_install_core_packages.R`:

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
pheatmap
ggpubr
cowplot
```

Heavy optional packages are checked by `packages.R` but not installed by default:

```text
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
```

SCENIC-related scripts also require the corresponding motif-ranking database files. These database files are not included and should be prepared following the SCENIC workflow.

## Output

The scripts generate `.rds`, `.pdf`, `.png`, and `.csv` files. Most outputs are written to `output/`; a few plotting scripts write files to the working directory because they follow the original figure-export code.

Downstream scripts often require objects generated by earlier scripts. If a script reports a missing Seurat object, regulon AUC matrix, or predicted drug-response object, run the upstream script first or place the required object in the expected path.

## Repository scope

This repository documents the single-cell computational analysis only. It does not include raw sequencing matrices, large intermediate RDS files, clinical datasets, IHC images, wet-lab experimental data, or manual figure assembly files.
