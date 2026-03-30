read_rds_input <- function(input_dir, assay = "RNA") {
  library(Matrix)
  library(data.table)
  
  rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(rds_files) == 0)
    stop("No .rds file found in ", input_dir)
  if (length(rds_files) > 1)
    stop("Multiple .rds files found: ", paste(rds_files, collapse = ", "))
  
  obj <- readRDS(rds_files[1])
  
  # Optional user-provided metadata
  metadata_file <- list.files(input_dir, pattern = "_metadata\\.(csv|tsv)$", full.names = TRUE)
  use_metadata <- length(metadata_file) == 1
  
  # -------------------------
  # SingleCellExperiment
  # -------------------------
  if (inherits(obj, "SingleCellExperiment")) {
    message("Detected SingleCellExperiment")
    
    counts <- SingleCellExperiment::counts(obj)
    counts <- as(counts, "dgCMatrix")
    
    # Ensure genes x cells
    if (nrow(counts) < ncol(counts)) {
      message("Counts appears transposed. Fixing orientation.")
      counts <- t(counts)
    }
    
    # Cells
    if (use_metadata) {
      cells <- fread(metadata_file[1])
      message("Using user-provided metadata: ", basename(metadata_file[1]))
    } else {
      cells <- as.data.table(SummarizedExperiment::colData(obj))
    }
    cells[, cell_id := colnames(counts)]
    
    # Genes
    genes <- as.data.table(SummarizedExperiment::rowData(obj))
    genes[, gene_id := rownames(counts)]
  }
  
  # -------------------------
  # Seurat
  # -------------------------
  else if (inherits(obj, "Seurat")) {
    message("Detected Seurat object")
    
    if (!assay %in% names(obj@assays))
      stop("Assay '", assay, "' not found in Seurat object")
    
    counts <- Seurat::GetAssayData(obj, assay = assay, layer = "counts")
    counts <- as(counts, "dgCMatrix")
    
    # Detect malformed Seurat object (cells as rows)
    if (all(grepl("^[ACGT]+$", rownames(counts)))) {
      message("Counts appears flipped (barcodes as rows). Fixing orientation.")
      counts <- t(counts)
    }
    
    # Cells
    if (use_metadata) {
      cells <- fread(metadata_file[1])
      message("Using user-provided metadata: ", basename(metadata_file[1]))
    } else {
      cells <- as.data.table(obj@meta.data)
    }
    cells[, cell_id := colnames(counts)]
    
    # Genes
    genes <- data.table(
      gene_id = rownames(counts)
    )
  }
  
  else {
    stop("Unsupported RDS object class: ", paste(class(obj), collapse = ", "))
  }
  
  # -------------------------
  # Final sanity checks
  # -------------------------
  if (ncol(counts) != nrow(cells))
    stop("Mismatch: ncol(counts) != nrow(cells)")
  
  if (nrow(counts) != nrow(genes))
    stop("Mismatch: nrow(counts) != nrow(genes)")
  
  if (anyDuplicated(colnames(counts)))
    stop("Duplicate cell names detected")
  
  if (anyDuplicated(rownames(counts)))
    stop("Duplicate gene names detected")
  
  list(
    counts = counts,   # genes x cells
    cells  = cells,
    genes  = genes
  )
}


