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

The runner skips optional heavy steps when their packages or input files are not available. This avoids a capsule failure caused by packages such as inferCNV, CopyKAT, SCENIC, or Monocle3.

For a full run after all dependencies and data have been prepared:

```bash
RUN_HEAVY_STEPS=true STRICT_CORE=true Rscript run_codeocean.R
```

The original strict entry point is still available:

```bash
Rscript run_all.R
```

## Data placement

Raw count matrices and large intermediate RDS files are not stored in the GitHub repository.

Recommended Data paths:

```text
/data/GSE138709/
/data/scRNA1_preprocessed.rds
/data/scRNA1_annotated.rds
/data/scRNA1_with_SCENIC_AUC.rds
/data/scRNA1_with_DrugPredictions.rds
```

Equivalent paths under `/data/output/` are also supported.

If no input data are found, `run_codeocean.R` exits without error and writes:

```text
output/CODEOCEAN_INPUT_REQUIRED.txt
```

This makes the missing-data reason explicit instead of producing a long dependency-related crash.

## Package strategy

Use `packages.R` as a package availability check:

```bash
Rscript packages.R
```

It writes:

```text
output/package_check.csv
```

System libraries can be installed with:

```bash
bash codeocean_system_deps.sh
```

The script installs system libraries only. R packages should be managed by the Code Ocean environment where possible.

## Practical submission setup

1. Import this GitHub repository into the Code section.
2. Upload processed RDS objects to the Data section.
3. Set the capsule run command to:

```bash
Rscript run_codeocean.R
```

4. Check `output/codeocean_run_log.csv` after running.
5. Use the original scripts and `run_all.R` as the full reproducibility record.

The stable runner is intended to make the capsule inspectable and executable even when heavyweight optional steps are not rerun inside the review environment.
