# Code Ocean notes

This repository contains the code used for the single-cell RNA-seq analyses. The Code Ocean capsule should contain both the code and the required input data or processed intermediate objects.

## Recommended command

Use the reviewer-safe runner as the main command:

```bash
Rscript run_codeocean.R
```

This command writes a log to:

```text
output/codeocean_run_log.csv
```

The runner uses the ordered script entry points `01_*.R` to `18_*.R`. It skips optional heavy steps when their packages or input files are not available, which prevents failures caused by inferCNV, CopyKAT, SCENIC, or Monocle3.

For a full run after all dependencies and data have been prepared:

```bash
RUN_HEAVY_STEPS=true STRICT_CORE=true Rscript run_codeocean.R
```

The original strict entry point is still available:

```bash
Rscript run_all.R
```

## Data placement

GSE138709 is provided by GEO as processed UMI count matrices in `GSE138709_RAW.tar`. The archive contains eight files ending with `_UMI.csv.gz`.

Recommended Data paths:

```text
/data/GSE138709_RAW.tar
/data/GSE138709_RAW/
/data/scRNA1_preprocessed.rds
/data/scRNA1_annotated.rds
/data/scRNA1_with_SCENIC_AUC.rds
/data/scRNA1_with_DrugPredictions.rds
```

Equivalent paths under `/data/output/` are also supported for processed RDS objects.

If `GSE138709_RAW.tar` is present, `run_codeocean.R` extracts it into:

```text
GSE138709_RAW/
```

The preprocessing script then reads the extracted UMI CSV files directly.

If no input data are found, `run_codeocean.R` exits without error and writes:

```text
output/CODEOCEAN_INPUT_REQUIRED.txt
```

This makes the missing-data reason explicit instead of producing a long dependency-related crash.

## Package strategy

The standard installation and run order is documented in `INSTALL.md`.

In a clean Code Ocean capsule, use:

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

`00_install_core_packages.R` installs only core CRAN packages used by the default reviewer-safe workflow. `packages.R` checks both core and optional package availability and writes:

```text
output/package_check.csv
```

Do not install Ubuntu `r-cran-*` or `r-bioc-*` packages in this capsule unless the base image is specifically designed for that.

## Practical submission setup

1. Import this GitHub repository into the Code section.
2. Upload `GSE138709_RAW.tar` or processed RDS objects to the Data section.
3. Set the capsule run command to:

```bash
Rscript run_codeocean.R
```

4. Check `output/codeocean_run_log.csv` after running.
5. Use the ordered entry scripts and `run_all.R` as the full reproducibility record.

The stable runner is intended to make the capsule inspectable and executable even when heavyweight optional steps are not rerun inside the review environment.
