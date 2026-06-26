# YB-1 ICC Single-Cell — Code Ocean Reproduction Package

Single-cell analysis of intrahepatic cholangiocarcinoma (iCCA, GSE138709), reproducing all
single-cell figures from **Figure 1 / Figure 3 / Figure S1 / Figure S3**. All outputs are
written to `results/` with article-figure-number prefixes (PDF + PNG; FigS2 is PDF only).

> Figure 2, Figure 4–7, Figure S4–S7, and Tables 1–2 correspond to clinical/IHC/survival/
> wet-lab/tabular data and are outside the scope of this package.

---

## 1. Directory Structure

```
YB_1_scRNA/
├── 00_load_data.R            # Source first: loads all objects + sets output directory
├── RUN_ORDER.md              # Execution order + which object each figure uses
├── Fig1A_cell_annotation.R
├── Fig1B … FigS3E_to_H.R    # 26 figure scripts (each self-contained)
├── data/                     # Processed objects (~4.4 GB); not tracked in git
└── results/                  # Output figures (PDF + PNG); generated at run time
```

## 2. Environment

R ≥ 4.2. Standard CRAN packages: `Seurat`, `ggplot2`, `dplyr`, `tidyr`, `reshape2`,
`aplot`, `scales`, `patchwork`, `colorspace`, `ggbeeswarm`, `ggrepel`, `RColorBrewer`,
`pheatmap`, `ks`, `igraph`.

**Bioconductor / special packages:**
- `Nebulosa` — kernel-density UMAPs (Fig1C, FigS3C, FigS3D, Fig3F)
- `ComplexHeatmap` + `circlize` — tile heatmaps (FigS3B, FigS3E–H)
- `infercnv` (FigS2) — no binary on current Bioc; must be **compiled from source**.
  Its dependency `argparse` requires Python at load time; set
  `PYTHON3=/path/to/python` before installing (Python needs the `argparse` and `json` modules).
- `monocle` **v2** (Fig1H–K pseudotime) — removed from current Bioc; install from the
  **Bioc 3.17 archive**: first install `HSMMSingleCell` (source) + `VGAM`, `DDRTree`,
  `fastICA`, then install `monocle_2.28.0.tar.gz`.

**monocle2 compatibility with newer dplyr / igraph:**
- If the environment pins **dplyr 1.0.x + igraph 2.0.x** (the monocle2-era versions),
  Fig1H–K runs without modification.
- With newer versions (dplyr ≥ 1.1, igraph ≥ 2.1), `Fig1H_to_1K_trajectory.R` applies
  **three compatibility patches automatically**:
  ① SE verbs (`group_by_` / `select_` / …) → modern equivalents;
  ② `igraph::dfs(neimode=)` → `mode=`;
  ③ `nei()` → `.nei()` inside monocle internals.
  Patches activate only when the newer versions are detected; they are skipped silently
  in older environments.

## 3. How to Run

```r
setwd("<repo-root>")          # or set env var YB1_DATA to the data folder
source("00_load_data.R")      # load once per session
source("Fig1A_cell_annotation.R")   # source individual scripts, or see RUN_ORDER.md
```

**Fig1H–K note:** the first monocle2 DDRTree run takes ~24 minutes (12,290 cells, high
peak memory). The script supports caching: set `YB1_CDS_CACHE=/some/path/cds.rds` and
the trajectory object is saved after the first run, making subsequent re-plots
instantaneous. Leave this unset for a full recompute in the Code Ocean capsule.

---

## 4. Figure Reproduction Status — All Complete, Verified Against Published Figures

| Figure | Content | Object | Script |
|--------|---------|--------|--------|
| 1A | Cell-type annotation UMAP | Full dataset | `Fig1A_cell_annotation.R` |
| 1B | Sample origin UMAP | Full dataset | `Fig1B_sample_origin_UMAP.R` |
| 1C | YB-1 density UMAP | Full dataset | `Fig1C_YBX1_density.R` |
| 1D | YB-1 expression by Tumor/Adjacent | Full dataset | `Fig1D_YBX1_by_sample.R` |
| 1E | CopyKAT UMAP | Epithelial subset | `Fig1E_copykat_UMAP.R` |
| 1F | YB-1 expression UMAP | Epithelial subset | `Fig1F_YBX1_expression_UMAP.R` |
| 1G | YB-1 by aneuploid/diploid | Epithelial subset | `Fig1G_YBX1_boxdensity.R` |
| 1H–K | Monocle2 pseudotime (Pseudotime / CopyKAT / State / YB-1) | Epithelial subset | `Fig1H_to_1K_trajectory.R` |
| 3A | SCENIC regulon heatmap (per-cell) | scRNAauc + regulonAUC | `Fig3A_SCENIC_regulon_heatmap.R` |
| 3B | YB-1 regulon AUC (UMAP + raincloud) | scRNAauc | `Fig3B_…` + `Fig1G_…` (bottom panels) |
| 3C | GSVA YB-1 target genes (UMAP + raincloud) | scRNAauc | `Fig3C_…` + `Fig1G_…` (bottom panels) |
| 3D | Cisplatin sensitive/resistant — triple raincloud | scRNAauc_cis | `Fig3D_cisplatin_boxdensity.R` |
| 3E | Gemcitabine — triple raincloud | scRNAauc + gemcitabine groups | `Fig3E_gemcitabine_boxdensity.R` |
| 3F | ABCB1 / ABCC1 / ABCC2 / MVP density UMAPs | Full dataset | `Fig3F_ABC_transporter_density.R` |
| S1A | Cluster UMAP | Full dataset | `FigS1A_cluster_UMAP.R` |
| S1B | Cell-type marker violin plots | Full dataset | `FigS1B_celltype_marker_violin.R` |
| S1C | Cell-type proportions | Full dataset | `FigS1C_celltype_proportion.R` |
| S1D | Epithelial sub-cluster UMAP | Epithelial subset | `FigS1D_cholangiocyte_subcluster_UMAP.R` |
| S1E | Subtype UMAP (LPC / Chol / Malig) | Epithelial subset | `FigS1E_subtype_UMAP.R` |
| S1F | Subtype marker violin plots | Epithelial subset | `FigS1F_subtype_marker_violin.R` |
| S2 | inferCNV heatmap | infercnv_obj | `FigS2_inferCNV_heatmap.R` |
| S3A | Regulon volcano plot (cisplatin) | scRNAauc_cis + regulonAUC | `FigS3A_regulon_volcano.R` |
| S3B | Regulon AUC 2-column heatmap (cisplatin) | scRNAauc_cis + regulonAUC | `FigS3B_regulon_heatmap.R` |
| S3C | Cisplatin lnIC50 density UMAP | scRNAauc_cis | `FigS3C_cisplatin_IC50_density.R` |
| S3D | YB-1 regulon AUC density UMAP | scRNAauc | `FigS3D_YBX1_regulon_density.R` |
| S3E–H | YB-1/ABC 5-gene × group tile heatmaps | Full / subset / cisplatin / gemcitabine | `FigS3E_to_H_ABC_heatmap.R` |

---

## 5. Data Manifest (`data/`, ~4.4 GB)

Large data files are not tracked in this repository. Place them in `data/` before running,
or set the environment variable `YB1_DATA` to their location.

| File | Contents |
|------|----------|
| `integrated_FindClusters_umap_singleR_0.5_curate_cnv_malignant_2.Rdata` | Full-dataset main object scRNA1 (33,990 cells) |
| `cholangiocyte_subset_FindClusters_sct_umap.Rdata` | Epithelial subset scRNA_sub (12,290 cells, own UMAP) |
| `scRNAauc.rds` | Subset + SCENIC regulon AUC (YBX1_extended_907g) + GSVA YB-1 target scores |
| `scRNAauc_with_CisplatinGroup_Cisplatin.rds` | scRNAauc + lnIC50 + cisplatin sensitivity groups |
| `3.4_regulonAUC.Rds` / `4.1_binaryRegulonActivity.Rds` / `2.5_regulonTargetsInfo.Rds` | SCENIC regulon matrices and target-gene table |
| `scRNA_aneuploid_with_Cisplatin_groups.rds` / `…Gemcitabine….rds` | Aneuploid subset + drug-sensitivity groups (Fig3D/3E, FigS3G/S3H) |
| `infercnv_obj.rds` | Final inferCNV object (FigS2) |
| `normal_copykat_clustering_results.rds` | CopyKAT raw clustering output |

`00_load_data.R` reads these files and exposes them as `scRNA1`, `scRNA_sub`, `scRNAauc`,
`scRNAauc_cis`, `scRNA_cis`, `scRNA_gem`, `regulonAUC_mat`, and `binaryReg_mat` for use
by the figure scripts.
