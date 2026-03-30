#!/usr/bin/env Rscript

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(edgeR)
  library(Matrix)
  library(rmarkdown)
  library(dplyr)
})

# ------------------------------------------------------------------
# Arguments
# ------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: 05_reporting.R <base_dir> <output_dir>")

base_dir    <- args[1]
reports_dir <- args[2]

dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------
# Find inputs
# ------------------------------------------------------------------
cellsets <- list.dirs(file.path(base_dir, "outputs", "cell_sets"), recursive = FALSE, full.names = FALSE)

# Load results_index
index_path <- file.path(base_dir, "outputs", "index", "results_index.tsv")
results_index <- if (file.exists(index_path)) fread(index_path) else data.table()

# ------------------------------------------------------------------
# Helper function: summarize a cell set
# ------------------------------------------------------------------
summarize_cellset <- function(cellset) {
  cat("Processing cell set:", cellset, "\n")
  
  # Make per-cell set directory
  cs_dir <- file.path(reports_dir, cellset)
  dir.create(cs_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Paths
  summary_path <- file.path(base_dir, "outputs", "cell_sets", cellset, "summary.tsv")
  qc_path      <- file.path(base_dir, "outputs","pseudobulk", cellset, "qc.tsv")
  
  # Read QC and summary
  qc      <- if (file.exists(qc_path)) fread(qc_path) else data.table()
  summary <- if (file.exists(summary_path)) fread(summary_path) else data.table()
  
  # Collect DE results
  cell_dge_dirs <- list.dirs(file.path(base_dir,"outputs","dge", cellset), recursive = FALSE, full.names = TRUE)
  de_results <- list()
  for (d in cell_dge_dirs) {
    result_file <- file.path(d, "result.tsv.gz")
    if (file.exists(result_file)) {
      de_results[[basename(d)]] <- fread(result_file)
    }
  }
  
  # Render per-cell set HTML
  report_file <- file.path(cs_dir, paste0(cellset, "_report.html"))
  rmarkdown::render(
    input = file.path(base_dir, "scripts", "templates", "cellset_report.Rmd"),
    output_file = report_file,
    params = list(
      cellset = cellset,
      summary = summary,
      qc = qc,
      de_results = de_results
    ),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )
  
  # Light TSV summary: top 20 genes per contrast
  summary_out <- data.table()
  for (contrast in names(de_results)) {
    top20 <- head(de_results[[contrast]][order(FDR)], 20)
    top20[, contrast := contrast]
    summary_out <- rbind(summary_out, top20, fill = TRUE)
  }
  fwrite(summary_out, file.path(cs_dir, paste0(cellset, "_top20.tsv")), sep = "\t")
}

# ------------------------------------------------------------------
# Loop over all cell sets
# ------------------------------------------------------------------
for (cs in cellsets) {
  summarize_cellset(cs)
}

# ------------------------------------------------------------------
# Global summary
# ------------------------------------------------------------------
global_de <- list()
for (cs in cellsets) {
  cell_dge_dirs <- list.dirs(file.path(base_dir, "outputs", "dge", cs), recursive = FALSE, full.names = TRUE)
  for (d in cell_dge_dirs) {
    result_file <- file.path(d, "result.tsv.gz")
    if (file.exists(result_file)) {
      de <- fread(result_file)
      de[, cell_set := cs]
      de[, contrast := basename(d)]
      global_de[[paste(cs, basename(d), sep = "_")]] <- de
    }
  }
}
global_de_dt <- rbindlist(global_de, fill = TRUE)

# Write global summary TSV
fwrite(global_de_dt, file.path(reports_dir, "global_de_summary.tsv"), sep = "\t")

# Optional global HTML report
global_report_file <- file.path(reports_dir, "global_report.html")
rmarkdown::render(
  input = file.path(base_dir, "scripts", "templates", "global_report.Rmd"),
  output_file = global_report_file,
  params = list(
    de = global_de_dt,
    results_index = results_index
  ),
  envir = new.env(parent = globalenv()),
  quiet = TRUE
)

cat("Reporting complete! HTML and TSV outputs written to:", reports_dir, "\n")
