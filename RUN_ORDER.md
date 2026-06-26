# Run Order

All figures are written to `../results/` with article-figure-number prefixes (PDF + PNG;
FigS2 is PDF only). For environment setup and dependencies (monocle2, infercnv, Fig1H–K
caching), see Section 2 of `../README.md`.

## 1. Load data first (once per session)

```r
source("00_load_data.R")
```

Loads: `scRNA1` (full dataset, 33,990 cells), `scRNA_sub` (epithelial subset, 12,290 cells),
`scRNAauc` (subset + regulon AUC / GSVA), `scRNAauc_cis` (+ cisplatin groups),
`scRNA_cis`, `scRNA_gem`, `regulonAUC_mat`, `binaryReg_mat`.

## 2. Run all 26 figures at once (order-independent; each script is self-contained)

```r
source("00_load_data.R")
for (s in c(
  "Fig1A_cell_annotation.R", "Fig1B_sample_origin_UMAP.R", "Fig1C_YBX1_density.R",
  "Fig1D_YBX1_by_sample.R", "Fig1E_copykat_UMAP.R", "Fig1F_YBX1_expression_UMAP.R",
  "Fig1G_YBX1_boxdensity.R", "Fig1H_to_1K_trajectory.R",
  "Fig3A_SCENIC_regulon_heatmap.R", "Fig3B_YBX1_regulon_AUC_UMAP.R", "Fig3C_GSVA_YBX1_targets_UMAP.R",
  "Fig3D_cisplatin_boxdensity.R", "Fig3E_gemcitabine_boxdensity.R", "Fig3F_ABC_transporter_density.R",
  "FigS1A_cluster_UMAP.R", "FigS1B_celltype_marker_violin.R", "FigS1C_celltype_proportion.R",
  "FigS1D_cholangiocyte_subcluster_UMAP.R", "FigS1E_subtype_UMAP.R", "FigS1F_subtype_marker_violin.R",
  "FigS2_inferCNV_heatmap.R",
  "FigS3A_regulon_volcano.R", "FigS3B_regulon_heatmap.R", "FigS3C_cisplatin_IC50_density.R",
  "FigS3D_YBX1_regulon_density.R", "FigS3E_to_H_ABC_heatmap.R"
)) source(s)
```

> Verified: all 26 scripts run to completion in a clean session from a single loader call (26/26 OK).

## 3. Script → Figure → Object reference

| Script | Figure | Object |
|--------|--------|--------|
| Fig1A_cell_annotation.R | Fig 1A | scRNA1 |
| Fig1B_sample_origin_UMAP.R | Fig 1B | scRNA1 |
| Fig1C_YBX1_density.R | Fig 1C | scRNA1 |
| Fig1D_YBX1_by_sample.R | Fig 1D | scRNA1 |
| Fig1E_copykat_UMAP.R | Fig 1E | scRNA_sub |
| Fig1F_YBX1_expression_UMAP.R | Fig 1F | scRNA_sub |
| Fig1G_YBX1_boxdensity.R | Fig 1G (+ Fig 3B/3C bottom panels) | scRNAauc |
| Fig1H_to_1K_trajectory.R | Fig 1H–K | scRNA_sub → Monocle2 |
| Fig3A_SCENIC_regulon_heatmap.R | Fig 3A | regulonAUC_mat + scRNAauc |
| Fig3B_YBX1_regulon_AUC_UMAP.R | Fig 3B (top) | scRNAauc |
| Fig3C_GSVA_YBX1_targets_UMAP.R | Fig 3C (top) | scRNAauc |
| Fig3D_cisplatin_boxdensity.R | Fig 3D | scRNAauc_cis |
| Fig3E_gemcitabine_boxdensity.R | Fig 3E | scRNAauc + gemcitabine groups |
| Fig3F_ABC_transporter_density.R | Fig 3F | scRNA1 |
| FigS1A_cluster_UMAP.R | Fig S1A | scRNA1 |
| FigS1B_celltype_marker_violin.R | Fig S1B | scRNA1 |
| FigS1C_celltype_proportion.R | Fig S1C | scRNA1 |
| FigS1D_cholangiocyte_subcluster_UMAP.R | Fig S1D | scRNA_sub |
| FigS1E_subtype_UMAP.R | Fig S1E | scRNA_sub |
| FigS1F_subtype_marker_violin.R | Fig S1F | scRNA_sub |
| FigS2_inferCNV_heatmap.R | Fig S2 | infercnv_obj |
| FigS3A_regulon_volcano.R | Fig S3A | scRNAauc_cis + regulonAUC_mat |
| FigS3B_regulon_heatmap.R | Fig S3B | scRNAauc_cis + regulonAUC_mat |
| FigS3C_cisplatin_IC50_density.R | Fig S3C | scRNAauc_cis |
| FigS3D_YBX1_regulon_density.R | Fig S3D | scRNAauc |
| FigS3E_to_H_ABC_heatmap.R | Fig S3E–H | scRNA1 / scRNA_sub / scRNA_cis / scRNA_gem |
