detect_input_type <- function(input_dir) {
  files <- list.files(input_dir)
  files_lower <- tolower(files)
  
  # ----------------------------
  # RDS
  # ----------------------------
  if (any(grepl("\\.rds$", files_lower))) return("rds")
  
  # ----------------------------
  # H5AD
  # ----------------------------
  if (any(grepl("\\.h5ad$", files_lower))) return("h5ad")
  
  # ----------------------------
  # ZARR
  # ----------------------------
  if (any(grepl("\\.zarr$", files_lower))) return("zarr")
  
  # ----------------------------
  # Parse / Splitpipe
  # Requires: cell_metadata.csv, count_matrix.mtx, all_genes.csv
  # Optional prefixes allowed
  # ----------------------------
  if (any(grepl("cell_metadata\\.csv$", files_lower)) &&
      any(grepl("count_matrix\\.mtx(\\.gz)?$", files_lower)) &&
      any(grepl("all_genes\\.csv$", files_lower))) {
    return("parse_splitpipe")
  }
  
  # ----------------------------
  # 10X / matrix option
  # Requires: matrix.mtx(.gz), barcodes.tsv, genes.tsv
  # Optional prefixes allowed
  # ----------------------------
  if (any(grepl("matrix\\.mtx(\\.gz)?$", files_lower)) &&
      any(grepl("barcodes\\.tsv(\\.gz)?$", files_lower)) &&
      any(grepl("genes\\.tsv(\\.gz)?$", files_lower))) {
    return("standard_matrix")
  }
  
  stop("Could not detect input type in ", input_dir)
}

