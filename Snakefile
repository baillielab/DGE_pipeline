# ----------------------------
# Snakefile
# ----------------------------

import os
from pathlib import Path

configfile: "config.yaml"

BASE_DIR = Path(config["input"]["base_dir"])
DATA_DIR = BASE_DIR / config["input"].get("data_dir", "data")

# List of samples and cell sets
SAMPLES = [p.name for p in DATA_DIR.iterdir() if p.is_dir()]
CELL_SETS = list(config["cell_sets"].keys())

# ----------------------------
# DGE test configuration
# ----------------------------

CONTRASTS = config.get("dge", {}).get("tests", [])

# Build valid (cellset, contrast) combinations respecting per-test restrictions
DGE_COMBOS = []
for test in CONTRASTS:
    test_name = test["name"]
    # Use restricted cell sets if provided, otherwise all
    valid_cellsets = test.get("cell_sets", CELL_SETS)
    for cs in valid_cellsets:
        DGE_COMBOS.append((cs, test_name))

# DGE result files for rule all and S05_reporting
dge_files = [
    str(BASE_DIR / "outputs" / "dge" / cs / contrast / "result.tsv.gz")
    for cs, contrast in DGE_COMBOS
]

# ----------------------------
# Rule all
# ----------------------------

rule all:
    input:
        # Merged files
        str(BASE_DIR / "outputs" / "merged" / "genes.tsv"),
        str(BASE_DIR / "outputs" / "merged" / "cells.parquet"),
        str(BASE_DIR / "outputs" / "merged" / "counts.mtx.gz"),

        # Cell set membership
        expand(
            str(BASE_DIR / "outputs" / "cell_sets" / "{cellset}" / "membership.parquet"),
            cellset=CELL_SETS
        ),

        # Pseudobulk
        expand(
            str(BASE_DIR / "outputs" / "pseudobulk" / "{cellset}" / "counts.mtx.gz"),
            cellset=CELL_SETS
        ),

        # DGE results (only allowed combos)
        *dge_files,

        # Per-cellset reports
        expand(
            str(BASE_DIR  / "reports" / "{cellset}" / "{cellset}_report.html"),
            cellset=CELL_SETS
        ),
        expand(
            str(BASE_DIR  / "reports" / "{cellset}" / "{cellset}_top20.tsv"),
            cellset=CELL_SETS
        ),

        # Global report
        str(BASE_DIR  / "reports" / "global_report.html"),
        str(BASE_DIR  / "reports" / "global_de_summary.tsv")

# ----------------------------
# S01 Ingest
# ----------------------------
rule S01_ingest:
    input:
        config="config.yaml"
    output:
        genes=str(BASE_DIR / "outputs" / "{sample}" / "S01_ingest" / "genes.tsv"),
        cells=str(BASE_DIR / "outputs" / "{sample}" / "S01_ingest" / "cells.parquet"),
        counts=str(BASE_DIR / "outputs" / "{sample}" / "S01_ingest" / "counts.mtx.gz")
    params:
        sample="{sample}"
    shell:
        "Rscript scripts/01_ingest.R {input.config} {params.sample}"


# ----------------------------
# S01 Merge
# ----------------------------
rule S01_merge:
    input:
        expand(str(BASE_DIR / "outputs" / "{sample}" / "S01_ingest" / "genes.tsv"), sample=SAMPLES),
        expand(str(BASE_DIR / "outputs" / "{sample}" / "S01_ingest" / "cells.parquet"), sample=SAMPLES),
        expand(str(BASE_DIR / "outputs" / "{sample}" / "S01_ingest" / "counts.mtx.gz"), sample=SAMPLES),
        config="config.yaml"
    output:
        genes=str(BASE_DIR / "outputs" / "merged" / "genes.tsv"),
        cells=str(BASE_DIR / "outputs" / "merged" / "cells.parquet"),
        counts=str(BASE_DIR / "outputs" / "merged" / "counts.mtx.gz")
    shell:
        "Rscript scripts/01.1_merge.R {input.config}"


# ----------------------------
# S02 Resolve Cell Sets
# ----------------------------
rule S02_resolve_cell_sets:
    input:
        config="config.yaml",
        cells=str(BASE_DIR / "outputs" / "merged" / "cells.parquet")
    output:
        membership=expand(
            str(BASE_DIR / "outputs" / "cell_sets" / "{cellset}" / "membership.parquet"),
            cellset=CELL_SETS
        ),
        summary=expand(
            str(BASE_DIR / "outputs" / "cell_sets" / "{cellset}" / "summary.tsv"),
            cellset=CELL_SETS
        )
    shell:
        "Rscript scripts/02_resolve_cell_sets.R {input.config}"

# ----------------------------
# S03 Build Pseudobulk
# ----------------------------
rule S03_build_pseudobulk:
    input:
        config="config.yaml",
        counts=str(BASE_DIR / "outputs" / "merged" / "counts.mtx.gz"),
        cells=str(BASE_DIR / "outputs" / "merged" / "cells.parquet"),
        genes=str(BASE_DIR / "outputs" / "merged" / "genes.tsv"),
        membership=str(BASE_DIR / "outputs" / "cell_sets" / "{cellset}" / "membership.parquet")
    output:
        counts=str(BASE_DIR / "outputs" / "pseudobulk" / "{cellset}" / "counts.mtx.gz"),
        units=str(BASE_DIR / "outputs" / "pseudobulk" / "{cellset}" / "units.tsv"),
        qc=str(BASE_DIR / "outputs" / "pseudobulk" / "{cellset}" / "qc.tsv"),
        genes=str(BASE_DIR / "outputs" / "pseudobulk" / "{cellset}" / "genes_used.tsv")
    shell:
        "Rscript scripts/03_pseudobulk.R {input.config} {wildcards.cellset}"


# ----------------------------
# S04 Run DGE
# ----------------------------
rule S04_run_dge:
    input:
        counts=str(BASE_DIR / "outputs" / "pseudobulk" / "{cellset}" / "counts.mtx.gz"),
        units=str(BASE_DIR / "outputs" / "pseudobulk" / "{cellset}" / "units.tsv"),
        genes_used=str(BASE_DIR / "outputs" / "pseudobulk" / "{cellset}" / "genes_used.tsv"),
        config="config.yaml"
    output:
        result=str(BASE_DIR / "outputs" / "dge" / "{cellset}" / "{contrast}" / "result.tsv.gz"),
        model_json=str(BASE_DIR / "outputs" / "dge" / "{cellset}" / "{contrast}" / "model.json"),
        diagnostics=str(BASE_DIR / "outputs" / "dge" / "{cellset}" / "{contrast}" / "diagnostics.tsv")
    params:
        cellset="{cellset}",
        contrast="{contrast}"
    shell:
        "Rscript scripts/04_DGE.R {input.config} {params.cellset} {params.contrast}"

# ----------------------------
# S05 Reporting
# ----------------------------
rule S05_reporting:
    input:
        # Per-cell set summaries
        expand(str(BASE_DIR / "outputs" / "cell_sets" / "{cellset}" / "summary.tsv"), cellset=CELL_SETS),
        expand(str(BASE_DIR / "outputs" / "pseudobulk" / "{cellset}" / "qc.tsv"), cellset=CELL_SETS),
        # Only allowed DGE results
        expand(str(BASE_DIR / "outputs" / "dge" / "{cellset}" / "{contrast}" / "result.tsv.gz"),
               zip,
               cellset=[cs for cs, _ in DGE_COMBOS],
               contrast=[c for _, c in DGE_COMBOS])
    output:
        # Each cell set gets its own folder
        reports=expand(str(BASE_DIR  / "reports" / "{cellset}" / "{cellset}_report.html"), cellset=CELL_SETS),
        top20=expand(str(BASE_DIR  / "reports" / "{cellset}" / "{cellset}_top20.tsv"), cellset=CELL_SETS),
        # Global report
        global_report=str(BASE_DIR  / "reports" / "global_report.html"),
        global_tsv=str(BASE_DIR  / "reports" / "global_de_summary.tsv")
    params:
        reports_dir=str(BASE_DIR  / "reports")
    shell:
        "Rscript scripts/05_report.R {BASE_DIR} {params.reports_dir}"
