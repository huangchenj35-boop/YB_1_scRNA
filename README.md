# YB-1 ICC single-cell — Code Ocean 复现包

GSE138709 肝内胆管癌(iCCA)单细胞分析,复现文章 **Figure 1 / Figure 3 / Figure S1 / Figure S3**
中的全部单细胞图。所有图输出到 `results/`,文件名带文章图号前缀(PDF + PNG;FigS2 为 PDF)。

> Figure 2、Figure 4–7、Figure S4–S7、Table 1–2 为临床/IHC/生存/湿实验/表格,无对应单细胞代码,不在本包范围。

---

## 1. 目录结构

```
Code_Ocean/
├── code/
│   ├── 00_load_data.R        # 先 source:载入所有对象 + 设输出目录
│   ├── RUN_ORDER.md          # 运行顺序 + 每图用哪个对象
│   └── Fig1A … FigS3E_to_H   # 26 个最终图脚本(每个自包含)
├── data/                     # 处理好的对象(~4.4 GB,< 5 GB 上限)
│   └── DATA_MANIFEST.md
├── results/                  # 运行后出图(PDF + PNG)
└── README.md
```

## 2. 运行环境(镜像里要装)

R ≥ 4.2。常规 CRAN 包:`Seurat`、`ggplot2`、`dplyr`、`tidyr`、`reshape2`、`aplot`、`scales`、
`patchwork`、`colorspace`、`ggbeeswarm`、`ggrepel`、`RColorBrewer`、`pheatmap`、`ks`、`igraph`。

**Bioconductor / 特殊包:**
- `Nebulosa`(密度图:Fig1C、FigS3C、FigS3D、Fig3F)
- `ComplexHeatmap` + `circlize`(方格热图:FigS3B、FigS3E–H)
- `infercnv`(FigS2)—— 当前 Bioc 无 binary,需**源码编译**;其依赖 `argparse` 加载时要找 Python,
  装前设环境变量 `PYTHON3=/path/to/python`(Python 需含 argparse、json 模块)
- `monocle` **v2**(Fig1H–K 拟时序)—— 已从当前 Bioc 下架,从 **Bioc 3.17 归档源码**装:
  先装 `HSMMSingleCell`(3.17 实验数据包源码)+ `VGAM`/`DDRTree`/`fastICA` 等,再装 `monocle_2.28.0.tar.gz`

**monocle2 与新版 dplyr/igraph 的兼容:**
- 若镜像 pin 了 **dplyr 1.0.x + igraph 2.0.x**(monocle2 同年代版本),Fig1H–K 直接可跑。
- 若用新版(dplyr ≥ 1.1、igraph ≥ 2.1),`Fig1H_to_1K_trajectory.R` **内置 3 个兼容补丁**会自动生效:
  ① dplyr SE 动词(`group_by_`/`select_`/…)→ 现代动词;② igraph `dfs(neimode=)` → `mode=`;
  ③ monocle 函数里的 `nei()` → `.nei()`。补丁只在检测到新版时启用,老版环境自动跳过。

## 3. 怎么跑

```r
setwd("code")
source("00_load_data.R")            # 一次会话载入一次
source("Fig1A_cell_annotation.R")   # 逐个 source,或见 RUN_ORDER.md 的一键循环
```

**Fig1H–K 注意**:首次跑 monocle2 DDRTree 约 **24 分钟**(12290 细胞,内存峰值高)。脚本支持缓存:
设环境变量 `YB1_CDS_CACHE=/some/path/cds.rds`,第一次跑会把建好的轨迹对象存下,之后**秒级重画**。
镜像里首次复现不设缓存即可(完整重算)。

---

## 4. 图复现状态 —— 全部完成、逐张与原图核对一致

| 图 | 内容 | 对象 | 脚本 |
|---|---|---|---|
| 1A | 细胞注释 UMAP | 全细胞 | `Fig1A_cell_annotation.R` |
| 1B | 样本来源 UMAP | 全细胞 | `Fig1B_sample_origin_UMAP.R` |
| 1C | YB-1 密度 UMAP | 全细胞 | `Fig1C_YBX1_density.R` |
| 1D | YB-1 by Tumor/Adjacent | 全细胞 | `Fig1D_YBX1_by_sample.R` |
| 1E | CopyKAT UMAP | 上皮子集 | `Fig1E_copykat_UMAP.R` |
| 1F | YB-1 表达 UMAP | 上皮子集 | `Fig1F_YBX1_expression_UMAP.R` |
| 1G | YB-1 by aneuploid/diploid | 上皮子集 | `Fig1G_YBX1_boxdensity.R` |
| 1H–K | Monocle2 拟时序(Pseudotime/CopyKAT/State/YB-1) | 上皮子集 | `Fig1H_to_1K_trajectory.R` |
| 3A | SCENIC regulon 热图(per-cell) | scRNAauc + regulonAUC | `Fig3A_SCENIC_regulon_heatmap.R` |
| 3B | YB-1 regulon AUC(UMAP + raincloud) | scRNAauc | `Fig3B_…` + `Fig1G_…`(下) |
| 3C | GSVA YB-1 靶基因(UMAP + raincloud) | scRNAauc | `Fig3C_…` + `Fig1G_…`(下) |
| 3D | 顺铂敏感/耐药 三联 raincloud | scRNAauc_cis | `Fig3D_cisplatin_boxdensity.R` |
| 3E | 吉西他滨 三联 raincloud | scRNAauc + 吉西他滨分组 | `Fig3E_gemcitabine_boxdensity.R` |
| 3F | ABCB1/ABCC1/ABCC2/MVP 密度 UMAP | 全细胞 | `Fig3F_ABC_transporter_density.R` |
| S1A | 聚类 UMAP | 全细胞 | `FigS1A_cluster_UMAP.R` |
| S1B | 细胞类型 marker 小提琴 | 全细胞 | `FigS1B_celltype_marker_violin.R` |
| S1C | 细胞比例 | 全细胞 | `FigS1C_celltype_proportion.R` |
| S1D | 子集聚类 UMAP | 上皮子集 | `FigS1D_cholangiocyte_subcluster_UMAP.R` |
| S1E | 亚型 UMAP(LPC/Chol/Malig) | 上皮子集 | `FigS1E_subtype_UMAP.R` |
| S1F | 亚型 marker 小提琴 | 上皮子集 | `FigS1F_subtype_marker_violin.R` |
| S2 | inferCNV 热图 | infercnv_obj | `FigS2_inferCNV_heatmap.R` |
| S3A | regulon 火山(顺铂) | scRNAauc_cis + regulonAUC | `FigS3A_regulon_volcano.R` |
| S3B | regulon AUC 2列热图(顺铂) | scRNAauc_cis + regulonAUC | `FigS3B_regulon_heatmap.R` |
| S3C | 顺铂 lnIC50 密度 UMAP | scRNAauc_cis | `FigS3C_cisplatin_IC50_density.R` |
| S3D | YB-1 regulon 密度 UMAP | scRNAauc | `FigS3D_YBX1_regulon_density.R` |
| S3E–H | YB-1/ABC 5基因×分组 方格热图 | 全细胞/子集/顺铂/吉西他滨 | `FigS3E_to_H_ABC_heatmap.R` |

---

## 5. 数据清单(`data/`,~4.4 GB)

| 文件 | 用途 |
|---|---|
| `integrated_FindClusters_umap_singleR_0.5_curate_cnv_malignant_2.Rdata` | 全细胞主对象 scRNA1(33990) |
| `cholangiocyte_subset_FindClusters_sct_umap.Rdata` | 上皮子集 scRNA_sub(12290,自身 umap) |
| `scRNAauc.rds` | 子集 + SCENIC regulon AUC(YBX1_extended_907g)+ GSVA YB-1 靶基因 |
| `scRNAauc_with_CisplatinGroup_Cisplatin.rds` | scRNAauc + lnIC50 + 顺铂分组 |
| `3.4_regulonAUC.Rds` / `4.1_binaryRegulonActivity.Rds` / `2.5_regulonTargetsInfo.Rds` | SCENIC regulon 矩阵/靶基因表 |
| `scRNA_aneuploid_with_Cisplatin_groups.rds` / `…Gemcitabine…` | 非整倍体 + 药敏分组(S3G/S3H、Fig3D/3E) |
| `infercnv_obj.rds` | inferCNV 最终对象(FigS2) |
| `normal_copykat_clustering_results.rds` | copyKAT 结果 |

`00_load_data.R` 把它们载入为 `scRNA1` / `scRNA_sub` / `scRNAauc` / `scRNAauc_cis` /
`scRNA_cis` / `scRNA_gem` / `regulonAUC_mat` 等,供各图脚本调用。
