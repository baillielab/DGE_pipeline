#!/usr/bin/env Rscript
rm(list = ls())

suppressPackageStartupMessages({
  library(yaml)
  library(SingleCellExperiment)
  library(scuttle)
  library(Matrix)
  library(edgeR)
  library(arrow)
  library(data.table)
})

# ----------------------------
# Arguments
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: 03_pseudobulk.R <config.yaml> <cellset_id>")

config_path <- args[1]
cellset_id  <- args[2]

cat("=== S03_build_pseudobulk ===\n")
cat("Cell set:", cellset_id, "\n")

config <- yaml::read_yaml(config_path)
base_dir <- config$input$base_dir

# ----------------------------
# Pseudobulk parameters
# ----------------------------
pb_cfg <- config$pseudobulk

unit_col       <- pb_cfg$aggregation$unit_col
condition_col  <- pb_cfg$aggregation$condition_col
batch_col      <- pb_cfg$aggregation$batch_col
donor_col      <- pb_cfg$aggregation$donor_col

min_cells_per_unit        <- pb_cfg$filtering$units$min_cells_per_unit
min_total_counts_per_unit <- pb_cfg$filtering$units$min_total_counts_per_unit
use_group_for_filter      <- pb_cfg$filtering$genes$use_group_for_filter

downsample_enabled <- pb_cfg$downsampling$enabled
downsample_max     <- pb_cfg$downsampling$max_cells_per_unit
downsample_seed    <- pb_cfg$downsampling$seed

set.seed(downsample_seed)

# ----------------------------
# Paths
# ----------------------------
merged_dir      <- file.path(base_dir, "outputs", "merged")
cellset_dir     <- file.path(base_dir, "outputs", "cell_sets", cellset_id)
out_dir         <- file.path(base_dir, "outputs", "pseudobulk", cellset_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts_path     <- file.path(merged_dir, "counts.mtx.gz")
cells_path      <- file.path(merged_dir, "cells.parquet")
genes_path      <- file.path(merged_dir, "genes.tsv")
membership_path <- file.path(cellset_dir, "membership.parquet")

# ----------------------------
# Load metadata
# ----------------------------
cat("Loading cells metadata...\n")
cells <- arrow::read_parquet(cells_path)

required_cols <- c("cell_id", unit_col, donor_col, condition_col, batch_col)
missing_cols <- setdiff(required_cols, colnames(cells))
if (length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse=", "))

# ----------------------------
# Load cell set membership
# ----------------------------
cat("Loading cell set membership...\n")
membership <- arrow::read_parquet(membership_path)
cells <- cells[cells$cell_id %in% membership$cell_id, ]
if (nrow(cells) == 0) stop("No cells found for cell set: ", cellset_id)

# ----------------------------
# Load genes
# ----------------------------
cat("Loading genes...\n")
genes <- fread(genes_path)

# ----------------------------
# Load counts
# ----------------------------
cat("Loading counts matrix...\n")
counts <- Matrix::readMM(gzfile(counts_path))
full_cells <- arrow::read_parquet(cells_path)
keep_idx <- which(full_cells$cell_id %in% membership$cell_id)
counts <- counts[, keep_idx]
cells <- full_cells[keep_idx, ]

# ----------------------------
# Optional downsampling
# ----------------------------
if (downsample_enabled) {
  cat("Applying downsampling...\n")
  keep_cells <- unlist(lapply(split(seq_len(nrow(cells)), cells[[unit_col]]), function(idx) {
    if (length(idx) > downsample_max) sample(idx, downsample_max) else idx
  }))
  cells <- cells[keep_cells, ]
  counts <- counts[, keep_cells]
}

# ----------------------------
# Build SCE
# ----------------------------
colData <- as.data.frame(cells[, .(
  sample = get(unit_col),
  donor = get(donor_col),
  condition = get(condition_col),
  batch = get(batch_col)
)])
colnames(colData) <- c(unit_col, donor_col, condition_col, batch_col)

sce <- SingleCellExperiment(
  assays = list(counts = as(counts, "dgCMatrix")),
  colData = colData
)

# ----------------------------
# Aggregate pseudobulk
# ----------------------------
cat("Aggregating counts per unit...\n")
agg <- scuttle::aggregateAcrossCells(sce, ids = sce[[unit_col]])
pb_counts <- assay(agg, "counts")
unit_df <- as.data.frame(colData(agg))
unit_df$library_size <- Matrix::colSums(pb_counts)
unit_df$n_cells <- as.numeric(table(cells[[unit_col]])[unit_df[[unit_col]]])

# ----------------------------
# Filter units
# ----------------------------
cat("Filtering units...\n")
keep_units <- unit_df$n_cells >= min_cells_per_unit &
  unit_df$library_size >= min_total_counts_per_unit
if (sum(keep_units) == 0) stop("All units filtered out for cell set: ", cellset_id)

pb_counts <- pb_counts[, keep_units]
unit_df   <- unit_df[keep_units, ]

# ----------------------------
# Gene filtering
# ----------------------------
cat("Filtering genes with edgeR::filterByExpr...\n")
dge <- DGEList(counts = pb_counts)
if (use_group_for_filter) {
  group <- factor(unit_df[[condition_col]])
  keep_genes <- filterByExpr(dge, group = group)
} else {
  keep_genes <- filterByExpr(dge)
}

pb_counts <- pb_counts[keep_genes, ]
genes_used <- genes$gene_id[keep_genes]

# ----------------------------
# QC table
# ----------------------------
qc_df <- unit_df[, c(unit_col, donor_col, condition_col, batch_col, "n_cells", "library_size")]

# ----------------------------
# Write outputs
# ----------------------------
cat("Writing outputs...\n")
tmp <- tempfile(fileext = ".mtx")
pb_counts <- as(pb_counts, "dgCMatrix")
Matrix::writeMM(pb_counts, tmp)
system2("gzip", c("-c", tmp), stdout = file.path(out_dir, "counts.mtx.gz"))
unlink(tmp)

fwrite(unit_df, file.path(out_dir, "units.tsv"), sep = "\t")
fwrite(qc_df, file.path(out_dir, "qc.tsv"), sep = "\t")
fwrite(data.frame(gene_id = genes_used), file.path(out_dir, "genes_used.tsv"), sep = "\t")

cat("S03 complete for:", cellset_id, "\n")

