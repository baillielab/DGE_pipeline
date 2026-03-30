# Pseudobulk DGE Pipeline (Snakemake)

This repository contains the Snakemake workflow and scripts for running the ODAP pseudobulk differential gene expression (DGE) pipeline.

## Overview

This pipeline performs **pseudobulk differential gene expression analysis** from single-cell or single-nucleus RNA-seq data.

It aggregates counts across cells to generate sample-level expression profiles for defined cell populations, then performs differential expression testing between conditions while accounting for donor effects and optional covariates.

The workflow is implemented in **Snakemake** for reproducibility and scalability on HPC systems.

## Pipeline Summary

The pipeline includes:

* Data ingestion from multiple formats
* Metadata harmonisation and mapping
* Cell set definition (metadata-driven)
* Pseudobulk aggregation (per sample/donor)
* Differential expression analysis (configurable model)
* Reporting (QC metrics, PCA, DE results)

## Input Support

Supported input formats include:

* Seurat (`.rds`)
* SingleCellExperiment (`.rds`)
* Parse Splitpipe outputs
* Standard Matrix Market (`.mtx`, `.tsv`)
* *(Planned)* AnnData (`.h5ad`)

**SEE [FULL DOCUMENTATION](https://git.ecdf.ed.ac.uk/odap-users-guide/odap-users-guide/-/wikis/Analysis-in-ODAP/Differential-Gene-Expression-Pipeline) FOR INPUT FORMAT REQUIREMENTS**
Multiple datasets and formats can be mixed within a single analysis.

## Key Concepts

### Cell Sets

Cells can be grouped into user-defined populations using:

* **Logical expressions** (e.g. metadata filters)
* **Set algebra** (e.g. intersections, unions)

### Pseudobulk

Counts are aggregated per:

* **sample (default)**
* optionally donor or other grouping variables

This enables robust statistical testing using bulk RNA-seq methods.

### Differential Expression

Flexible model design allows:

* specification of reference condition
* inclusion of covariates (e.g. donor, batch)
* multiple contrasts and tests

## Usage

This pipeline is designed for HPC execution using Slurm.

Basic steps:

Edit configuration:

```bash id="6n5k1r"
nano config.yaml
```

Run pipeline:

```bash id="k5p1jz"
sbatch run_pipeline.slurm
```

## Configuration

All parameters are defined in `config.yaml`, including:

* input data and metadata mapping
* cell set definitions
* pseudobulk settings
* differential expression model and contrasts

The pipeline is highly flexible and supports heterogeneous input datasets.

## Outputs

* Intermediate files (`outputs/`)
* Differential expression results (`outputs/dge/`)

  * per cell set and contrast
  * `results.tsv.gz` for downstream analysis
* HTML reports (`results/`)

  * per cell set
  * global summaries
* Top gene summaries

## Documentation

Full user guide (including setup, metadata structure, and configuration details) is available in the ODAP documentation:

[Full documentation and usage instructions](https://git.ecdf.ed.ac.uk/odap-users-guide/odap-users-guide/-/wikis/Analysis-in-ODAP/Differential-Gene-Expression-Pipeline)

## Notes

* Designed for multi-sample, multi-condition studies
* Supports complex experimental designs
* Combines single-cell resolution with bulk statistical robustness

## Requirements

These are all available in the ODAP DGE Pipeline container.


* Snakemake
* Python
* Slurm (for HPC execution)
* R
  * Matrix
  * data.table
  * arrow
  * SingleCellExperiment
  * zellkonverter
  * yaml
  * rlang
  * jsonlite
  * Seurat
  * fs
  * SeuratObject
  * scuttle
  * edgeR
  * r.utils
  * ggplot2
  * rmarkdown
  * dplyr
  * UpsetR
  * tibble
  * readr
  * knitr
  * scales
  * ggrepel

## Author

**Kathryn Campbell**
on behalf of The ODAP Team

30th March 2026

