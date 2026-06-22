# Installation and run order

This project uses a reviewer-safe Code Ocean workflow. The default run does not require all heavy single-cell packages to be installed. Heavy steps are skipped unless their packages and input objects are available.

## Recommended Code Ocean order

Run the following commands from the repository root.

### 1. Optional system libraries

Run this only in a fresh Code Ocean capsule when system libraries are missing:

```bash
bash codeocean_system_deps.sh
```

This script installs system libraries only. It does not install Ubuntu `r-cran-*` or `r-bioc-*` packages.

### 2. Core R packages

Install the core CRAN packages needed by the default workflow:

```bash
Rscript 00_install_core_packages.R
```

This installs only the core packages used by the default reviewer-safe workflow:

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

The installation check is written to:

```text
output/core_package_install_check.csv
```

### 3. Package check

Check both core and optional package availability:

```bash
Rscript packages.R
```

The package check is written to:

```text
output/package_check.csv
```

### 4. Reviewer-safe analysis run

Run the Code Ocean stable runner:

```bash
Rscript run_codeocean.R
```

The run log is written to:

```text
output/codeocean_run_log.csv
```

## Standard one-block command

In a clean Code Ocean capsule, the standard command block is:

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

## Heavy optional packages

The following packages are not installed by the default installer:

```text
monocle3
infercnv
copykat
SCENIC
AUCell
RcisTarget
GENIE3
GSVA
ComplexHeatmap
Nebulosa
```

These packages are used by trajectory, CNV, SCENIC, GSVA, and regulon-density steps. In the default reviewer-safe workflow, scripts requiring missing heavy packages are skipped and recorded in `output/codeocean_run_log.csv`.

To run all heavy steps after installing all required packages and preparing the corresponding input objects, use:

```bash
RUN_HEAVY_STEPS=true STRICT_CORE=true Rscript run_codeocean.R
```

## Do not use these commands in the Code Ocean capsule

Avoid mixing Ubuntu R binary packages with the CRAN R runtime used by the capsule:

```bash
apt-get install r-cran-*
apt-get install r-bioc-*
```

Also avoid changing apt repositories inside the running capsule unless the base image explicitly requires it.
