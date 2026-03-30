#!/usr/bin/env Rscript

# ==================================================================
# S04_run_dge.R
# Differential expression per cell set using edgeR
# Fully spec-compliant result table and diagnostics
# ==================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(Matrix)
  library(edgeR)
  library(jsonlite)
})

# ------------------------------------------------------------------
# Arguments
# ------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: S04_run_dge.R <config.yaml> <cellset_id>")

config_path <- args[1]
cellset_id  <- args[2]

cat("=== S04_run_dge ===\n")
cat("Cell set:", cellset_id, "\n")

# ------------------------------------------------------------------
# Load configuration
# ------------------------------------------------------------------
config <- yaml::read_yaml(config_path)
base_dir <- config$input$base_dir

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
pb_dir      <- file.path(base_dir, "outputs", "pseudobulk", cellset_id)
dge_dir     <- file.path(base_dir, "outputs", "dge", cellset_id)
index_path  <- file.path(base_dir, "outputs", "index", "results_index.tsv")

dir.create(dge_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(index_path), recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(pb_dir, "counts.mtx.gz")
units_path  <- file.path(pb_dir, "units.tsv")
genes_path  <- file.path(pb_dir, "genes_used.tsv")

# ------------------------------------------------------------------
# Load pseudobulk
# ------------------------------------------------------------------
cat("Loading pseudobulk data...\n")
counts <- readMM(gzfile(counts_path))
counts <- as(counts, "dgCMatrix")
units  <- fread(units_path)
genes  <- fread(genes_path)

if (ncol(counts) != nrow(units)) stop("Mismatch between counts columns and units rows")

# ------------------------------------------------------------------
# Constraints
# ------------------------------------------------------------------
constraints <- config$dge$constraints
if (nrow(units) < constraints$min_total_units) {
  stop("Total units below min_total_units constraint")
}

# ------------------------------------------------------------------
# Apply reference levels and factorize
# ------------------------------------------------------------------
design_cfg <- config$dge$design
design_formula <- as.formula(design_cfg$formula)

if (!is.null(design_cfg$reference_levels)) {
  for (var in names(design_cfg$reference_levels)) {
    if (!(var %in% colnames(units))) stop(paste("Reference variable not found:", var))
    units[[var]] <- relevel(factor(units[[var]]), ref = design_cfg$reference_levels[[var]])
  }
}

for (v in all.vars(design_formula)) {
  if (v %in% colnames(units) && !is.numeric(units[[v]])) {
    units[[v]] <- factor(units[[v]])
  }
}

# ------------------------------------------------------------------
# Build design matrix and edgeR model
# ------------------------------------------------------------------
design <- model.matrix(design_formula, data = units)
cat("Design matrix columns:\n")
print(colnames(design))

dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

# ------------------------------------------------------------------
# Prepare results index
# ------------------------------------------------------------------
if (!file.exists(index_path)) {
  fwrite(data.frame(
    cell_set  = character(),
    test_name = character(),
    status    = character(),
    reason    = character()
  ), index_path, sep = "\t")
}
results_index <- fread(index_path)

# ------------------------------------------------------------------
# Loop through tests
# ------------------------------------------------------------------
for (test_def in config$dge$tests) {
  
  test_name <- test_def$name
  test_type <- test_def$type
  
  if (!is.null(test_def$cell_sets) && !(cellset_id %in% test_def$cell_sets)) {
    cat("Skipping test", test_name, "for cell set", cellset_id, "\n")
    next
  }
  
  cat("Running test:", test_name, "\n")
  test_dir <- file.path(dge_dir, test_name)
  dir.create(test_dir, recursive = TRUE, showWarnings = FALSE)
  
  tryCatch({
    # --------------------------
    # Contrast test
    # --------------------------
    if (test_type == "contrast") {
      contrast_vec <- makeContrasts(contrasts = test_def$contrast, levels = design)
      qlf <- glmQLFTest(fit, contrast = contrast_vec)
      group1_name <- gsub("condition", "", test_def$contrast)
      group2_name <- design_cfg$reference_levels$condition
      group1_idx <- which(units$condition == group1_name)
      group2_idx <- which(units$condition == group2_name)
    }
    
    # --------------------------
    # Multi-coefficient test
    # --------------------------
    else if (test_type == "coef_set") {
      coef_idx <- match(test_def$coefficients, colnames(design))
      if (any(is.na(coef_idx))) stop("One or more coefficients not found in design matrix")
      qlf <- glmQLFTest(fit, coef = coef_idx)
      group1_idx <- 1:nrow(units)
      group2_idx <- integer(0)
      group1_name <- paste(test_def$coefficients, collapse="_")
      group2_name <- NA
    }
    
    else stop("Unknown test type")
    
    # --------------------------
    # Extract top results
    # --------------------------
    top <- topTags(qlf, n = Inf)$table
    top$FDR <- p.adjust(top$PValue, method = constraints$multiple_testing)
    
    logcounts <- cpm(dge, log=TRUE)
    
    mean_expr_group1 <- rowMeans(logcounts[, group1_idx, drop=FALSE])
    mean_expr_group2 <- if(length(group2_idx) > 0) rowMeans(logcounts[, group2_idx, drop=FALSE]) else NA
    prop_units_expr_group1 <- rowMeans(counts[, group1_idx, drop=FALSE] > 0)
    prop_units_expr_group2 <- if(length(group2_idx) > 0) rowMeans(counts[, group2_idx, drop=FALSE] > 0) else NA
    n_units_group1 <- length(group1_idx)
    n_units_group2 <- length(group2_idx)
    
    result_table <- data.frame(
      gene_id = genes$gene_id,
      gene_symbol = if("gene_symbol" %in% colnames(genes)) genes$gene_symbol else NA,
      logFC   = if ("logFC" %in% colnames(top)) top$logFC else NA,
      PValue  = top$PValue,
      FDR     = top$FDR,
      logCPM  = top$logCPM,
      stat    = if ("F" %in% colnames(top)) top$F else top$LR,
      AveExpr_or_logCPM = top$logCPM,
      n_units_total = nrow(units),
      mean_expr_group1 = mean_expr_group1,
      mean_expr_group2 = mean_expr_group2,
      prop_units_expr_group1 = prop_units_expr_group1,
      prop_units_expr_group2 = prop_units_expr_group2,
      n_units_group1 = n_units_group1,
      n_units_group2 = n_units_group2
    )
    
    fwrite(result_table, file.path(test_dir, "result.tsv.gz"), sep = "\t")
    
    # --------------------------
    # Save model metadata
    # --------------------------
    write_json(list(
      formula   = design_cfg$formula,
      test_name = test_name,
      test_type = test_type,
      engine    = "edgeR_qlf"
    ), file.path(test_dir, "model.json"), auto_unbox = TRUE, pretty = TRUE)
    
    # --------------------------
    # Save diagnostics
    # --------------------------
    diagnostics <- data.frame(
      n_genes = nrow(counts),
      n_units = nrow(units),
      common_dispersion = dge$common.dispersion
    )
    fwrite(diagnostics, file.path(test_dir, "diagnostics.tsv"), sep = "\t")
    
    # --------------------------
    # Update results index
    # --------------------------
    results_index <- rbind(
      results_index,
      data.frame(
        cell_set  = cellset_id,
        test_name = test_name,
        status    = "SUCCESS",
        reason    = ""
      ),
      fill = TRUE
    )
    
    cat("SUCCESS:", test_name, "\n")
    
  }, error = function(e) {
    results_index <- rbind(
      results_index,
      data.frame(
        cell_set  = cellset_id,
        test_name = test_name,
        status    = "FAILED",
        reason    = e$message
      ),
      fill = TRUE
    )
    cat("FAILED:", test_name, "\n")
    cat("Reason:", e$message, "\n")
  })
}

# ------------------------------------------------------------------
# Save results index
# ------------------------------------------------------------------
fwrite(results_index, index_path, sep = "\t")
cat("S04 complete for:", cellset_id, "\n")
