# 运行顺序 / Run order

所有图输出到 `../results/`,文件名带文章图号前缀(PDF + PNG;FigS2 为 PDF)。
环境与依赖见 `../README.md` 第 2 节(尤其 monocle2 / infercnv 的装法与 Fig1H–K 缓存)。

## 1. 先加载数据(每次会话一次)

```r
source("00_load_data.R")
```

载入对象:`scRNA1`(全细胞 33990)、`scRNA_sub`(上皮子集 12290)、`scRNAauc`(子集+regulon AUC/GSVA)、
`scRNAauc_cis`(+顺铂分组)、`scRNA_gsva`、`scRNA_cis`、`scRNA_gem`、`regulonAUC_mat`、`binaryReg_mat`。

## 2. 一键跑全部 26 张(顺序无关,逐个独立)

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

> 已验证:以上 26 个脚本在干净会话里从 loader 一次性全跑通(26/26 OK)。

## 3. 脚本 → 图 → 对象 对照

| 脚本 | 图 | 对象 |
|---|---|---|
| Fig1A_cell_annotation.R | 图1A | scRNA1 |
| Fig1B_sample_origin_UMAP.R | 图1B | scRNA1 |
| Fig1C_YBX1_density.R | 图1C | scRNA1 |
| Fig1D_YBX1_by_sample.R | 图1D | scRNA1 |
| Fig1E_copykat_UMAP.R | 图1E | scRNA_sub |
| Fig1F_YBX1_expression_UMAP.R | 图1F | scRNA_sub |
| Fig1G_YBX1_boxdensity.R | 图1G(+3B/3C 下) | scRNAauc |
| Fig1H_to_1K_trajectory.R | 图1H–K | scRNA_sub → Monocle2 |
| Fig3A_SCENIC_regulon_heatmap.R | 图3A | regulonAUC_mat + scRNAauc |
| Fig3B_YBX1_regulon_AUC_UMAP.R | 图3B 上 | scRNAauc |
| Fig3C_GSVA_YBX1_targets_UMAP.R | 图3C 上 | scRNAauc |
| Fig3D_cisplatin_boxdensity.R | 图3D | scRNAauc_cis |
| Fig3E_gemcitabine_boxdensity.R | 图3E | scRNAauc + 吉西他滨分组 |
| Fig3F_ABC_transporter_density.R | 图3F | scRNA1 |
| FigS1A_cluster_UMAP.R | 图S1A | scRNA1 |
| FigS1B_celltype_marker_violin.R | 图S1B | scRNA1 |
| FigS1C_celltype_proportion.R | 图S1C | scRNA1 |
| FigS1D_cholangiocyte_subcluster_UMAP.R | 图S1D | scRNA_sub |
| FigS1E_subtype_UMAP.R | 图S1E | scRNA_sub |
| FigS1F_subtype_marker_violin.R | 图S1F | scRNA_sub |
| FigS2_inferCNV_heatmap.R | 图S2 | infercnv_obj |
| FigS3A_regulon_volcano.R | 图S3A | scRNAauc_cis + regulonAUC_mat |
| FigS3B_regulon_heatmap.R | 图S3B | scRNAauc_cis + regulonAUC_mat |
| FigS3C_cisplatin_IC50_density.R | 图S3C | scRNAauc_cis |
| FigS3D_YBX1_regulon_density.R | 图S3D | scRNAauc |
| FigS3E_to_H_ABC_heatmap.R | 图S3E–H | scRNA1 / scRNA_sub / scRNA_cis / scRNA_gem |
