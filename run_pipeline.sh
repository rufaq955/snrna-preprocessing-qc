#!/usr/bin/env bash
#
# run_pipeline.sh - run all preprocessing and QC steps in order.
#
# Stops on the first failing step so a broken stage does not silently
# propagate empty or partial output into the next one.
#
# Usage:
#   bash run_pipeline.sh --samples config/samples.tsv --outdir data [--cuda]
#
set -euo pipefail

SAMPLES="config/samples.tsv"
OUTDIR="data"
WORKERS=4
CUDA_FLAG=""
SKIP_CONVERT=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samples) SAMPLES="$2"; shift 2 ;;
    --outdir)  OUTDIR="$2"; shift 2 ;;
    --workers) WORKERS="$2"; shift 2 ;;
    --cuda)    CUDA_FLAG="--cuda"; shift ;;
    --skip-convert) SKIP_CONVERT=1; shift ;;
    -h|--help) sed -n '2,11p' "$0"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

[[ -f "$SAMPLES" ]] || { echo "ERROR: sample manifest not found: $SAMPLES" >&2; exit 1; }

mkdir -p "$OUTDIR"/{interim/{h5ad,cellbender,r_h5,doublets,meta},qc,logs}

banner() { printf '\n\033[1m== %s ==\033[0m\n' "$1"; }

if [[ $SKIP_CONVERT -eq 0 ]]; then
  banner "Step 01 - 10x to AnnData"
  python scripts/01_10x_to_h5ad.py \
    --samples "$SAMPLES" \
    --outdir "$OUTDIR/interim/h5ad" \
    --workers "$WORKERS"
else
  banner "Step 01 - skipped (--skip-convert)"
fi

banner "Step 02 - CellBender ambient RNA removal"
python scripts/02_cellbender_remove_background.py \
  --indir "$OUTDIR/interim/h5ad" \
  --outdir "$OUTDIR/interim/cellbender" \
  $CUDA_FLAG

banner "Step 03 - Repack HDF5 for R"
bash scripts/03_ptrepack_for_r.sh \
  --indir "$OUTDIR/interim/cellbender" \
  --outdir "$OUTDIR/interim/r_h5"

banner "Step 04 - DoubletFinder"
Rscript scripts/04_doubletfinder.R \
  --indir "$OUTDIR/interim/r_h5" \
  --outdir "$OUTDIR/interim/doublets"

banner "Step 05a - QC metrics before filtering"
Rscript scripts/05_qc_filter.R \
  --indir "$OUTDIR/interim/doublets" \
  --outdir "$OUTDIR/qc" \
  --mode before

banner "Step 05b - QC filtering"
Rscript scripts/05_qc_filter.R \
  --indir "$OUTDIR/interim/doublets" \
  --outdir "$OUTDIR/qc" \
  --mode after

banner "Pipeline complete"
echo "Filtered objects : $OUTDIR/qc/*_qc.rds"
echo "QC summaries     : $OUTDIR/qc/qc_summary_*.csv"
echo "Doublet summary  : $OUTDIR/interim/doublets/doublet_summary.csv"
echo "Logs             : $OUTDIR/logs/"
