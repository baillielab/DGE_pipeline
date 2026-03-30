#!/usr/bin/env Rscript
rm(list = ls())

# ----------------------------
# Packages
# ----------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(Matrix)
  library(yaml)
  library(SeuratObject)
  library(Seurat)
})

# ----------------------------
# Arguments
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: 01_ingest.R <config.yaml> <sample_name>")
}

config <- yaml::read_yaml(args[1])
sample_name <- args[2]

# ----------------------------
# Config & paths
# ----------------------------
base <- config$input$base_dir
if (is.null(base)) stop("config$input$base_dir is required")

ingest_cfg <- config$ingest
if (is.null(ingest_cfg$default)) stop("ingest$default section is required")

# Input folder for this sample
input_dir <- file.path(base, "data", sample_name)
if (!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir)

# Pick mapping: per-sample override or default
mapping <- ingest_cfg$samples[[sample_name]]
if (is.null(mapping)) mapping <- ingest_cfg$default

required_fields <- c("sample_id_col", "donor_id_col", "condition_col",
                     "batch_col", "raw_cell_id_cols")
missing <- setdiff(required_fields, names(mapping))
if (length(missing)) stop("Missing ingest mapping fields: ", paste(missing, collapse=", "))

# Output dir
out_dir <- file.path(base, "outputs", sample_name, "S01_ingest")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Base directory: ", base)
message("Input directory: ", input_dir)
message("Sample name: ", sample_name)
message("Output directory: ", out_dir)

# ----------------------------
# Source helper functions
# ----------------------------
rfunc_dir <- file.path(base, "scripts", "R_functions")
source(file.path(rfunc_dir, "detect_input_type.R"))
source(file.path(rfunc_dir, "read_parse_inputs.R"))
source(file.path(rfunc_dir, "read_10x_inputs.R"))
source(file.path(rfunc_dir, "read_zarr_inputs.R"))
source(file.path(rfunc_dir, "read_rds_inputs.R"))
source(file.path(rfunc_dir, "read_h5ad_inputs.R"))

# ----------------------------
# Harmonisation
# ----------------------------
harmonise_cells <- function(cells, mapping, raw_cell_id_col) {
  setDT(cells)
  
  required <- c(mapping$sample_id_col,
                mapping$donor_id_col,
                mapping$condition_col,
                mapping$batch_col)
  missing <- setdiff(required, colnames(cells))
  if (length(missing)) stop("Missing metadata columns: ", paste(missing, collapse=", "))
  
  # Resolve raw cell IDs
  if (raw_cell_id_col == "rownames") {
    cells[, .raw_cell_id := rownames(cells)]
  } else {
    if (!raw_cell_id_col %in% colnames(cells)) {
      stop("Raw cell ID column not found: ", raw_cell_id_col)
    }
    cells[, .raw_cell_id := get(raw_cell_id_col)]
  }
  
  # Construct cell_id
  cells[, cell_id := paste(get(mapping$sample_id_col), .raw_cell_id, sep="_")]
  
  # Harmonise standard columns
  cells[, `:=`(
    sample_id = get(mapping$sample_id_col),
    donor_id  = get(mapping$donor_id_col),
    condition = get(mapping$condition_col),
    batch     = get(mapping$batch_col)
  )]
  
  if (anyDuplicated(cells$cell_id)) stop("cell_id is not unique after harmonisation")
  
  cells[, .raw_cell_id := NULL]
  cells[]
}

# ----------------------------
# Validation
# ----------------------------
validate_inputs <- function(counts, genes, cells) {
  if (!inherits(counts, "dgCMatrix")) stop("Counts must be a dgCMatrix")
  if (ncol(counts) != nrow(cells)) stop("Counts columns != number of cells")
  if (nrow(counts) != nrow(genes)) stop("Counts rows != number of genes")
  if (anyDuplicated(cells$cell_id)) stop("cell_id not unique")
  invisible(TRUE)
}

# ----------------------------
# Detect & read input
# ----------------------------
input_type <- detect_input_type(input_dir)
message("Detected input type: ", input_type)

dat <- switch(
  input_type,
  parse_splitpipe = read_parse_inputs(input_dir),
  standard_matrix = read_10x_inputs(input_dir),
  rds             = read_rds_input(input_dir),
  h5ad            = read_h5ad_input(input_dir),
  zarr            = read_zarr_input(input_dir),
  stop("Unsupported input type: ", input_type)
)

counts <- dat$counts
cells  <- dat$cells
genes  <- dat$genes

# ----------------------------
# Harmonise
# ----------------------------
raw_cell_id_col <- mapping$raw_cell_id_cols[[input_type]]
if (is.null(raw_cell_id_col)) stop("No raw_cell_id_col defined for input type: ", input_type)

cells <- harmonise_cells(cells, mapping, raw_cell_id_col)

rownames(counts) <- genes$gene_id
colnames(counts) <- cells$cell_id

validate_inputs(counts, genes, cells)

# ----------------------------
# Derived metrics
# ----------------------------
cells[, n_umis := Matrix::colSums(counts)[cell_id]]

# ----------------------------
# Outputs
# ----------------------------
fwrite(genes, file.path(out_dir, "genes.tsv"), sep="\t")
write_parquet(cells, file.path(out_dir, "cells.parquet"))

tmp <- tempfile(fileext = ".mtx")
Matrix::writeMM(counts, tmp)
system2("gzip", c("-c", tmp), stdout = file.path(out_dir, "counts.mtx.gz"))
unlink(tmp)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("S01_ingest completed successfully for sample: ", sample_name)
