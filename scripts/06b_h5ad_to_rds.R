#!/usr/bin/env Rscript
#
# Step 06b - Convert .h5ad to a Seurat .rds object (with metadata).
#
# Three-stage hand-off:
#   1. Convert()          h5ad -> h5seurat
#   2. LoadH5Seurat()     h5seurat -> Seurat object (counts only)
#   3. re-attach metadata from the CSV produced by step 06a
#
# The cell-name identity check before attaching is deliberate: a silent
# reordering here would mislabel every cell in the object.
#
# Usage:
#   Rscript scripts/06b_h5ad_to_rds.R --h5ad input.h5ad \
#     --metadata data/interim/meta/input_metadata.csv --out input.rds
#
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(data.table)
  library(optparse)
})

# `%||%` is base R from 4.4.0 onward; define it for older interpreters.
if (!exists("%||%", envir = baseenv())) `%||%` <- function(x, y) if (is.null(x)) y else x

option_list <- list(
  make_option("--h5ad", type = "character", help = "Input .h5ad file"),
  make_option("--metadata", type = "character", default = NULL,
              help = "Metadata CSV from step 06a (optional)"),
  make_option("--out", type = "character", default = NULL,
              help = "Output .rds path [default: <input>.rds]")
)
opts <- parse_args(OptionParser(option_list = option_list))

if (is.null(opts$h5ad)) stop("--h5ad is required")
if (!file.exists(opts$h5ad)) stop(sprintf("Not found: %s", opts$h5ad))

out_path <- opts$out %||% sub("\\.h5ad$", ".rds", opts$h5ad)
h5seurat_path <- tempfile(fileext = ".h5seurat")

message("Converting ", opts$h5ad, " -> h5seurat")
Convert(opts$h5ad, dest = h5seurat_path, overwrite = TRUE, verbose = FALSE)

message("Loading as Seurat object")
seurat_obj <- LoadH5Seurat(h5seurat_path, misc = FALSE, meta.data = FALSE)

if (!is.null(opts$metadata)) {
  message("Attaching metadata from ", opts$metadata)
  meta <- as.data.frame(fread(opts$metadata))
  rownames(meta) <- meta[[1]]
  meta[[1]] <- NULL

  if (identical(colnames(seurat_obj), rownames(meta))) {
    seurat_obj@meta.data <- meta
    message("  metadata attached (", ncol(meta), " columns)")
  } else if (setequal(colnames(seurat_obj), rownames(meta))) {
    meta <- meta[colnames(seurat_obj), , drop = FALSE]
    seurat_obj@meta.data <- meta
    message("  metadata attached after reordering to match cell order")
  } else {
    warning("Cell names do not match between data and metadata; ",
            "saving object WITHOUT metadata. Check step 06a output.")
  }
}

saveRDS(seurat_obj, file = out_path)
unlink(h5seurat_path)

message("Saved ", out_path, " (", ncol(seurat_obj), " cells x ",
        nrow(seurat_obj), " genes)")
