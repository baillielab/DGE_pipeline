read_zarr_input <- function(input_dir) {
  library(Matrix)
  library(data.table)
  library(zellkonverter)
  library(SingleCellExperiment)

  zarr_dirs <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)
  zarr_dirs <- zarr_dirs[grepl("\\.zarr$", zarr_dirs)]
  if (length(zarr_dirs) == 0) stop("No .zarr directory found in ", input_dir)
  if (length(zarr_dirs) > 1) stop("Multiple .zarr directories found: ", paste(zarr_dirs, collapse=", "))

  sce <- readZarr(zarr_dirs[1])
  if (!inherits(sce, "SingleCellExperiment")) stop("Failed to convert Zarr to SCE")

  counts <- counts(sce)
  counts <- as(counts, "dgCMatrix")
  counts <- t(counts)
  if (nrow(counts) < ncol(counts)) counts <- t(counts)

  cells <- as.data.table(colData(sce))
  if (!"cell_id" %in% colnames(cells)) {
    cells[, cell_id := rownames(colData(sce))]
  }

  genes <- as.data.table(rowData(sce))
  if (!"gene_id" %in% colnames(genes)) genes[, gene_id := rownames(rowData(sce))]

  list(counts = counts, cells = cells, genes = genes)
}
