rm(list = ls())

# ----------------------------
# Load packages
# ----------------------------
suppressPackageStartupMessages({
  library(yaml)
  library(arrow)
  library(data.table)
  library(rlang)
  library(jsonlite)
})

# ----------------------------
# Arguments
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: 02_resolve_cell_sets.R <config.yaml>")
}

# ----------------------------
# File paths
# ----------------------------
config <- yaml::read_yaml(args[1])
base <- config[["input"]][["base_dir"]]
input_dir <- file.path(base, "outputs", "merged")
output_base <- file.path(base, "outputs", "cell_sets")
dir.create(output_base, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Load metadata
# ----------------------------
cells <- as.data.table(read_parquet(file.path(input_dir, "cells.parquet")))

# Required columns
required_cols <- c("cell_id", "condition", "n_umis", "donor_id")
missing <- setdiff(required_cols, names(cells))
if (length(missing) > 0) {
  stop("cells.parquet is missing required columns: ", paste(missing, collapse = ", "))
}

# Extract cell set definitions
cell_sets_cfg <- config$cell_sets
cell_set_names <- names(cell_sets_cfg)

# ----------------------------
# Functions
# ----------------------------

## Safe evaluation of 'where' expressions
safe_where_eval <- function(expr_string, dt) {
  expr <- rlang::parse_expr(expr_string)
  allowed_ops <- c(">", "<", ">=", "<=", "==", "!=", "&", "|", "!", "%in%")
  
  has_illegal_call <- function(x) {
    if (is.call(x)) {
      fname <- as.character(x[[1]])
      if (!(fname %in% allowed_ops)) return(TRUE)
    }
    if (is.pairlist(x) || is.expression(x)) {
      return(any(vapply(x, has_illegal_call, logical(1))))
    }
    FALSE
  }
  
  if (has_illegal_call(expr)) stop("Only comparison and boolean operators are allowed in where expressions")
  res <- rlang::eval_tidy(expr, data = dt)
  if (!is.logical(res) || length(res) != nrow(dt)) stop("where expression must return a logical vector")
  which(res)
}

## Tokenize algebra expression
tokenize_algebra <- function(expr_string) {
  tokens <- unlist(strsplit(expr_string, "\\s+"))
  operators <- c("AND", "OR", "SETDIFF")
  list(
    values = tokens[!(tokens %in% operators)],
    ops    = tokens[tokens %in% operators]
  )
}

## Evaluate algebra left-to-right
eval_algebra <- function(values, ops, resolved_sets) {
  if (length(values) < 2) stop("Algebra must have at least two cell sets")
  result <- resolved_sets[[values[1]]]
  for (i in seq_along(ops)) {
    rhs <- resolved_sets[[values[i + 1]]]
    op  <- ops[i]
    if (is.null(rhs)) stop("Undefined cell set in algebra: ", values[i + 1])
    result <- switch(
      op,
      "AND"     = intersect(result, rhs),
      "OR"      = union(result, rhs),
      "SETDIFF" = setdiff(result, rhs),
      stop("Unsupported algebra operator: ", op)
    )
  }
  sort(result)
}

# ----------------------------
# Resolve cell sets
# ----------------------------
resolved_sets <- list()

# 1. Where-based sets
for (id in cell_set_names) {
  cfg <- cell_sets_cfg[[id]]
  if (!is.null(cfg$where)) {
    resolved_sets[[id]] <- sort(safe_where_eval(cfg$where, cells))
  }
}

# 2. Algebra-based sets
for (id in cell_set_names) {
  cfg <- cell_sets_cfg[[id]]
  if (!is.null(cfg$algebra)) {
    tokens <- tokenize_algebra(cfg$algebra)
    resolved_sets[[id]] <- eval_algebra(tokens$values, tokens$ops, resolved_sets)
  }
}

# ----------------------------
# Write outputs
# ----------------------------
for (id in names(resolved_sets)) {
  idx <- resolved_sets[[id]]
  out_dir <- file.path(output_base, id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Membership
  write_parquet(data.table(cell_id = cells$cell_id[idx]),
                file.path(out_dir, "membership.parquet"))
  
  # Definition
  write_json(cell_sets_cfg[[id]],
             file.path(out_dir, "definition.json"),
             pretty = TRUE, auto_unbox = TRUE)
  
  # Summary (per condition)
  sub <- cells[idx]
  summary <- sub[, .(
    n_cells  = .N,
    n_umis   = sum(n_umis, na.rm = TRUE),
    n_donors = uniqueN(donor_id)
  ), by = condition]
  
  fwrite(summary,
         file.path(out_dir, "summary.tsv"),
         sep = "\t")
  
  # Session info
  sink(file.path(out_dir, "sessionInfo.txt"))
  sessionInfo()
  sink()
}


