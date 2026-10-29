#!/usr/bin/env Rscript
#
# Step 04 - Doublet detection and removal with DoubletFinder.
#
# For each sample: normalise -> variable features -> scale -> PCA -> estimate
# the expected doublet count from the recovered-cell number -> sweep pK ->
# classify -> keep singlets only.
#
# The doublet rate (default 0.075, i.e. 7.5%) is the 10x Chromium expectation
# for roughly 5,000 recovered cells; override with --doublet-rate.
#
# Usage:
#   Rscript scripts/04_doubletfinder.R \
#     --indir data/interim/r_h5 --outdir data/interim/doublets \
#     --doublet-rate 0.075 --pcs 15
#
suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(optparse)
})

option_list <- list(
  make_option("--indir", type = "character", default = "data/interim/r_h5"),
  make_option("--outdir", type = "character", default = "data/interim/doublets"),
  make_option("--doublet-rate", type = "double", default = 0.075,
              help = "Expected doublet fraction [default %default]"),
  make_option("--pcs", type = "integer", default = 15,
              help = "Number of PCs for the pK sweep [default %default]"),
  make_option("--pn", type = "double", default = 0.25,
              help = "DoubletFinder pN [default %default]")
)
opts <- parse_args(OptionParser(option_list = option_list))

dir.create(opts$outdir, showWarnings = FALSE, recursive = TRUE)

files <- Sys.glob(file.path(opts$indir, "*.h5"))
if (length(files) == 0) {
  stop(sprintf("No .h5 files found in %s", opts$indir))
}

message(sprintf("Processing %d sample(s) from %s", length(files), opts$indir))

for (file_path in files) {
  base_name <- tools::file_path_sans_ext(basename(file_path))
  message(sprintf("\n=== %s ===", base_name))

  result <- tryCatch({
    counts <- Read10X_h5(file_path)
    obj <- CreateSeuratObject(counts = counts, project = base_name)

    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = opts$pcs, verbose = FALSE)

    # Expected number of doublets scales with recovered cell count.
    n_cells <- ncol(obj)
    n_exp <- round(opts$doublet_rate * n_cells)

    # pK is data-dependent: sweep and take the pK maximising BCmetric.
    sweep_list <- paramSweep(obj, PCs = 1:opts$pcs, sct = FALSE)
    sweep_stats <- summarizeSweep(sweep_list, GT = FALSE)
    bcmvn <- find.pK(sweep_stats)
    pk_best <- as.numeric(as.character(
      bcmvn$pK[which.max(bcmvn$BCmetric)]
    ))

    obj <- doubletFinder(
      obj,
      PCs = 1:opts$pcs,
      pN = opts$pn,
      pK = pk_best,
      nExp = n_exp,
      reuse.pANN = NULL,
      sct = FALSE
    )

    classification_col <- grep("DF.classifications", colnames(obj@meta.data),
                               value = TRUE)
    if (length(classification_col) == 0) {
      stop("DoubletFinder classification column not found")
    }

    n_doublets <- sum(obj@meta.data[[classification_col[1]]] == "Doublet")
    singlet_obj <- obj[, obj@meta.data[[classification_col[1]]] == "Singlet"]

    outfile <- file.path(opts$outdir, paste0(base_name, "_Doublet.rds"))
    saveRDS(singlet_obj, file = outfile)

    message(sprintf(
      "  cells: %d -> %d  |  pK: %.4f  |  doublets removed: %d  |  saved: %s",
      n_cells, ncol(singlet_obj), pk_best, n_doublets, basename(outfile)
    ))

    data.frame(
      sample = base_name,
      cells_before = n_cells,
      doublets_removed = n_doublets,
      cells_after = ncol(singlet_obj),
      pk = pk_best,
      expected_rate = opts$doublet_rate
    )
  }, error = function(e) {
    message(sprintf("  ERROR: %s", conditionMessage(e)))
    NULL
  })

  if (!is.null(result)) {
    summary_file <- file.path(opts$outdir, "doublet_summary.csv")
    write.table(result, summary_file,
                sep = ",", row.names = FALSE, col.names = !file.exists(summary_file),
                append = file.exists(summary_file), quote = FALSE)
  }
}

message("\nDone. Clean singlet objects written to ", opts$outdir)
