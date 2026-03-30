read_h5ad_input <- function(input_dir) {
  library(Matrix)
  library(data.table)
  library(zellkonverter)  # readH5AD
  library(SingleCellExperiment)
  
  h5ad_files <- list.files(input_dir, pattern = "\\.h5ad$", full.names = TRUE)
  if (length(h5ad_files) == 0) stop("No .h5ad file found in ", input_dir)
  if (length(h5ad_files) > 1) stop("Multiple .h5ad files found: ", paste(h5ad_files, collapse=", "))
  
  sce <- readH5AD(h5ad_files[1])
  if (!inherits(sce, "SingleCellExperiment")) stop("Failed to convert H5AD to SCE")
  
  assays_available <- assayNames(sce)
  
  if ("counts" %in% assays_available) {
    counts <- assay(sce, "counts")
  } else if ("X" %in% assays_available) {
    message("Using assay 'X' as counts")
    counts <- assay(sce, "X")
  } else {
    stop("No suitable assay found in H5AD file. Found: ",
         paste(assays_available, collapse = ", "))
  }
  
  counts <- as(counts, "dgCMatrix")
  
  # Ensure genes × cells orientation
  if (nrow(counts) < ncol(counts)) counts <- t(counts)
  
  # Optional user-provided metadata
  metadata_file <- list.files(input_dir, pattern = "_metadata\\.(csv|tsv)$", full.names = TRUE)
  use_metadata <- length(metadata_file) > 0
  
  if (use_metadata) {
    cells <- fread(metadata_file[1])
    message("Using user-provided metadata: ", basename(metadata_file[1]))
    if (!"cell_id" %in% colnames(cells)) {
      cells[, cell_id := colnames(counts)]
    }
  } else {
    cells <- as.data.table(colData(sce))
    if (!"cell_id" %in% colnames(cells)) {
      cells[, cell_id := colnames(sce)]
    }
  }
  
  genes <- as.data.table(rowData(sce))
  if (!"gene_id" %in% colnames(genes)) genes[, gene_id := rownames(rowData(sce))]
  
  list(
    counts = counts,
    cells  = cells,
    genes  = genes
  )
}
