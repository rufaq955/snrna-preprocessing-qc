#!/usr/bin/env bash
#
# Step 03 - Repack CellBender HDF5 output for R/Seurat.
#
# CellBender writes a filtered .h5 whose layout Seurat's Read10X_h5 does not
# always accept. Copying just the /matrix group through h5repack (ptrepack)
# with compression yields a compact file Seurat reads without complaint.
#
# Usage:
#   scripts/03_ptrepack_for_r.sh [--indir DIR] [--outdir DIR] [--complevel N]
#
set -euo pipefail

INDIR="data/interim/cellbender"
OUTDIR="data/interim/r_h5"
COMPLEVEL=5
SUFFIX="_cellbender_filtered.h5"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --indir)     INDIR="$2"; shift 2 ;;
    --outdir)    OUTDIR="$2"; shift 2 ;;
    --complevel) COMPLEVEL="$2"; shift 2 ;;
    -h|--help)   sed -n '2,12p' "$0"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

command -v ptrepack >/dev/null 2>&1 || {
  echo "ERROR: ptrepack not found. Install it with: pip install tables" >&2
  exit 1
}

mkdir -p "$OUTDIR"

shopt -s nullglob
files=("$INDIR"/*"$SUFFIX")

if [[ ${#files[@]} -eq 0 ]]; then
  echo "ERROR: no *${SUFFIX} files in ${INDIR}" >&2
  exit 1
fi

echo "Repacking ${#files[@]} file(s) -> ${OUTDIR} (complevel=${COMPLEVEL})"

for file in "${files[@]}"; do
  base="$(basename "$file" "$SUFFIX")"
  target="$OUTDIR/${base}_R.h5"

  if [[ -f "$target" ]]; then
    echo "  skip (exists): $(basename "$target")"
    continue
  fi

  echo "  repack: $(basename "$file") -> $(basename "$target")"
  ptrepack --complevel "$COMPLEVEL" "${file}:/matrix" "${target}:/matrix"
done

echo "Done. R-ready files in ${OUTDIR}"
