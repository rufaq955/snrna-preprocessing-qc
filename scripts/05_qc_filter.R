#!/usr/bin/env Rscript
#
# Step 05 - Cell-level QC filtering and mitochondrial/ribosomal annotation.
#
# Computes per-cell QC metrics and applies hard thresholds. Two modes:
#
#   --mode before   build objects and report metrics only (no filtering)
#   --mode after    apply the thresholds and save the filtered object
#
# Running both modes lets you produce the standard before/after QC figure
# that reviewers expect in a methods section.
#
# Default thresholds: nFeature_RNA > 200, nCount_RNA > 500,
#                     nCount_RNA < 10000, percent.mt <= 10
#
# Usage:
#   Rscript scripts/05_qc_filter.R --indir data/interim/doublets \
#     --outdir data/qc --mode after
#
suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
})

option_list <- list(
  make_option("--indir", type = "character", default = "data/interim/doublets"),
  make_option("--outdir", type = "character", default = "data/qc"),
  make_option("--mode", type = "character", default = "after",
              help = "'before' or 'after' [default %default]"),
  make_option("--min-features", type = "integer", default = 200),
  make_option("--min-counts", type = "integer", default = 500),
  make_option("--max-counts", type = "integer", default = 10000),
  make_option("--max-mt", type = "double", default = 10,
              help = "Maximum percent mitochondrial reads [default %default]")
)
opts <- parse_args(OptionParser(option_list = option_list))

if (!opts$mode %in% c("before", "after")) {
  stop("--mode must be either 'before' or 'after'")
}

dir.create(opts$outdir, showWarnings = FALSE, recursive = TRUE)

files <- Sys.glob(file.path(opts$indir, "*.rds"))
if (length(files) == 0) {
  stop(sprintf("No .rds files found in %s", opts$indir))
}

message(sprintf("Mode: %s | %d sample(s)", opts$mode, length(files)))

qc_rows <- list()

for (file_path in files) {
  # Sample names may carry a trailing suffix from earlier steps; strip it so
  # downstream labels stay consistent across the pipeline.
  sample_name <- gsub("_.*", "", tools::file_path_sans_ext(basename(file_path)))
  message(sprintf("  %s", sample_name))

  result <- tryCatch({
    data <- readRDS(file_path)
    counts <- GetAssayData(data, layer = "counts")
    obj <- CreateSeuratObject(counts = counts, project = sample_name)

    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")

    n_before <- ncol(obj)

    if (identical(opts$mode, "after")) {
      obj <- subset(
        obj,
        subset = nFeature_RNA > opts$min_features &
          nCount_RNA > opts$min_counts &
          nCount_RNA < opts$max_counts &
          percent.mt <= opts$max_mt
      )
      outfile <- file.path(opts$outdir, paste0(sample_name, "_qc.rds"))
      saveRDS(obj, file = outfile)
    } else {
      outfile <- NA_character_
    }

    message(sprintf(
      "    cells: %d -> %d | median genes: %.0f | median UMIs: %.0f | median %%mt: %.1f",
      n_before, ncol(obj),
      median(obj$nFeature_RNA), median(obj$nCount_RNA), median(obj$percent.mt)
    ))

    data.frame(
      sample = sample_name,
      mode = opts$mode,
      cells_before = n_before,
      cells_after = ncol(obj),
      median_genes = median(obj$nFeature_RNA),
      median_umis = median(obj$nCount_RNA),
      median_percent_mt = median(obj$percent.mt),
      output = outfile,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message(sprintf("    ERROR: %s", conditionMessage(e)))
    NULL
  })

  if (!is.null(result)) qc_rows[[length(qc_rows) + 1]] <- result
}

if (length(qc_rows) > 0) {
  summary_df <- do.call(rbind, qc_rows)
  summary_file <- file.path(opts$outdir, sprintf("qc_summary_%s.csv", opts$mode))
  write.csv(summary_df, summary_file, row.names = FALSE)
  message("\nSummary written to ", summary_file)
} else {
  message("\nNo samples processed successfully.")
}

message("Done.")
