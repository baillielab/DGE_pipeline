read_parse_inputs <- function(input_dir) {
  library(Matrix)
  library(data.table)
  
  # ----------------------------
  # Detect count matrix
  # ----------------------------
  mtx_files <- list.files(
    input_dir,
    pattern = "count.*\\.mtx(\\.gz)?$",
    full.names = TRUE
  )
  
  if (length(mtx_files) == 0) {
    stop("No count_matrix.mtx(.gz) file found in ", input_dir)
  }
  if (length(mtx_files) > 1) {
    stop("Multiple count_matrix files found in ", input_dir)
  }
  
  counts <- readMM(mtx_files[1])
  counts <- as(counts, "dgCMatrix")
  counts <- t(counts)  # genes x cells
  
  # ----------------------------
  # Detect genes
  # ----------------------------
  gene_files <- list.files(
    input_dir,
    pattern = "all_genes.*\\.csv$",
    full.names = TRUE
  )
  
  if (length(gene_files) == 0) {
    stop("No all_genes.csv file found in ", input_dir)
  }
  if (length(gene_files) > 1) {
    stop("Multiple all_genes files found in ", input_dir)
  }
  
  genes <- fread(gene_files[1])
  
  # Optional: standardise gene_id column if needed
  if (!"gene_id" %in% colnames(genes)) {
    setnames(genes, 1, "gene_id")
  }
  
  # ----------------------------
  # Detect metadata
  # ----------------------------
  metadata_file <- list.files(
    input_dir,
    pattern = "_metadata\\.(csv|tsv)$",
    full.names = TRUE
  )
  
  if (length(metadata_file) == 1) {
    cells <- fread(metadata_file[1])
    message("Using user-provided metadata: ", basename(metadata_file[1]))
    
  } else {
    # fallback: detect cell_metadata file
    cell_files <- list.files(
      input_dir,
      pattern = "cell_metadata.*\\.csv$",
      full.names = TRUE
    )
    
    if (length(cell_files) == 0) {
      stop("No cell_metadata.csv and no metadata file found in ", input_dir)
    }
    if (length(cell_files) > 1) {
      stop("Multiple cell_metadata files found in ", input_dir)
    }
    
    cells <- fread(cell_files[1])
    message("No user metadata found, using: ", basename(cell_files[1]))
  }
  
  # ----------------------------
  # Validate Parse requirement
  # ----------------------------
  if (!"bc_wells" %in% colnames(cells)) {
    stop("Parse cell metadata must contain 'bc_wells'")
  }
  
  list(
    counts = counts,
    genes  = genes,
    cells  = cells
  )
}

