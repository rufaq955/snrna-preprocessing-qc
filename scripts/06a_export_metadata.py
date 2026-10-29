#!/usr/bin/env python3
"""
Step 06a - Export cell-level metadata from .h5ad for the R hand-off.

SeuratDisk's `LoadH5Seurat(meta.data = FALSE)` is used in step 06b because
converting the full object in one pass is memory-hungry and sometimes drops
factor levels. Exporting `adata.obs` to CSV and re-attaching it in R keeps
the metadata intact and verifiable.

Usage
-----
    python scripts/06a_export_metadata.py --h5ad input.h5ad --outdir data/interim/meta
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=Path("data/interim/meta"))
    args = parser.parse_args()

    if not args.h5ad.exists():
        raise SystemExit(f"Not found: {args.h5ad}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.h5ad, backed="r")

    outfile = args.outdir / f"{args.h5ad.stem}_metadata.csv"
    adata.obs.to_csv(outfile)

    print(f"Wrote {outfile} ({adata.n_obs} cells x {adata.obs.shape[1]} columns)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
