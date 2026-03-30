read_10x_inputs <- function(input_dir) {
  library(Matrix)
  library(data.table)
  
  # ----------------------------
  # Detect count matrix
  # ----------------------------
  mtx_files <- list.files(
    input_dir,
    pattern = "matrix.*\\.mtx(\\.gz)?$",
    full.names = TRUE
  )
  
  if (length(mtx_files) == 0) {
    stop("No matrix.mtx(.gz) file found in ", input_dir)
  }
  if (length(mtx_files) > 1) {
    stop("Multiple matrix files found in ", input_dir)
  }
  
  counts <- readMM(mtx_files[1])
  counts <- as(counts, "dgCMatrix")

  # ----------------------------
  # Detect genes / features
  # ----------------------------
  genes_files <- list.files(
    input_dir,
    pattern = "(genes|features).*\\.tsv(\\.gz)?$",
    full.names = TRUE
  )
  
  if (length(genes_files) == 0) {
    stop("No genes/features file found in ", input_dir)
  }
  if (length(genes_files) > 1) {
    stop("Multiple genes/features files found in ", input_dir)
  }
  
  genes <- fread(genes_files[1], header = F)
  
  # Standardise to gene_id + gene_name
  if (ncol(genes) >= 2) {
    setnames(genes, 1:2, c("gene_id", "gene_name"))
  } else if (ncol(genes) == 1) {
    setnames(genes, 1, c("gene_id"))
  } else {
    stop("Genes file can have 1 or 2 columns only; gene_id and optional gene_name")
  }
  
  # ----------------------------
  # Detect optional metadata
  # ----------------------------
  metadata_file <- list.files(
    input_dir,
    pattern = "metadata\\.(csv|tsv)$",
    full.names = TRUE
  )
  
  if (length(metadata_file) == 1) {
    cells <- fread(metadata_file[1])
    message("Using user-provided metadata: ", basename(metadata_file[1]))
    
  } else {
    # fallback to barcodes
    barcode_files <- list.files(
      input_dir,
      pattern = "barcodes.*\\.tsv(\\.gz)?$",
      full.names = TRUE
    )
    
    if (length(barcode_files) == 0) {
      stop("No barcodes file found and no metadata provided")
    }
    if (length(barcode_files) > 1) {
      stop("Multiple barcode files found in ", input_dir)
    }
    
    cells <- fread(barcode_files[1], header = F)
    setnames(cells, 1, "barcode")
    
    message("No metadata file found, using barcodes only")
  }
  
  list(
    counts = counts,
    genes  = genes,
    cells  = cells
  )
}

