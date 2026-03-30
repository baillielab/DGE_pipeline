#!/usr/bin/env Rscript
rm(list = ls())

library(data.table)
library(arrow)
library(Matrix)
library(fs)

# ----------------------------
# Arguments
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)


# ----------------------------
# Config
# ----------------------------
library(yaml)

config <- yaml::read_yaml(args[1])
ingest_cfg <- config$ingest
base<-config[["input"]][["base_dir"]]

required_cols <- c("cell_id", "donor_id", "condition", "batch")

# ----------------------------
# Sample directories
# ----------------------------
data_dir <- file.path(base, "data")
samples <- dir_ls(data_dir, type = "directory")

out_dir <- file.path(base, "outputs", "merged")
dir_create(out_dir, recurse = TRUE)

all_cells <- list()
all_genes <- list()
all_counts <- list()

# ----------------------------
# Read S01 outputs per sample
# ----------------------------
for (sample in samples) {
  sample_name <- path_file(sample)
  s01_dir <- file.path(base, "outputs", sample_name, "S01_ingest")
  
  if (!dir_exists(s01_dir)) {
    warning("S01_ingest folder not found for sample: ", sample_name)
    next
  }
  
  # Read outputs
  cells <- read_parquet(file.path(s01_dir, "cells.parquet"))
  genes <- fread(file.path(s01_dir, "genes.tsv"))
  counts <- Matrix::readMM(file.path(s01_dir, "counts.mtx.gz"))
  
  # Prefix cell_id with sample to ensure uniqueness
  cells$cell_id <- paste0(sample_name, "_", cells$cell_id)
  colnames(counts) <- cells$cell_id
  
  # Check for missing required columns
  missing_cols <- setdiff(required_cols, colnames(cells))
  if (length(missing_cols)) {
    stop("Sample ", sample_name, " is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  
  # Track all extra columns
  extra_cols <- setdiff(colnames(cells), required_cols)
  if (length(extra_cols)) {
    message("Sample ", sample_name, " has extra metadata columns: ",
            paste(extra_cols, collapse = ", "))
  }
  
  all_cells[[sample_name]] <- as.data.table(cells)
  all_genes[[sample_name]] <- genes
  all_counts[[sample_name]] <- counts
}

# ----------------------------
# Merge outputs
# ----------------------------
# Cells
merged_cells <- rbindlist(all_cells, use.names = TRUE, fill = TRUE)

# Genes
# Assume all samples have same genes; take first
merged_genes <- all_genes[[1]]

# Counts
merged_counts <- do.call(cbind, all_counts)

# ----------------------------
# Validation
# ----------------------------
if (!all(required_cols %in% colnames(merged_cells))) {
  stop("Merged cells table is missing required columns")
}

if (ncol(merged_counts) != nrow(merged_cells)) {
  stop("Merged counts columns do not match number of merged cells")
}

if (nrow(merged_counts) != nrow(merged_genes)) {
  stop("Merged counts rows do not match number of genes")
}

# ----------------------------
# Save merged outputs
# ----------------------------
fwrite(merged_genes, file.path(out_dir, "genes.tsv"), sep = "\t")
write_parquet(merged_cells, file.path(out_dir, "cells.parquet"))

tmp <- tempfile(fileext = ".mtx")
Matrix::writeMM(merged_counts, tmp)
system2("gzip", c("-c", tmp), stdout = file.path(out_dir, "counts.mtx.gz"))
unlink(tmp)

message("Merged S01 outputs saved to: ", out_dir)

