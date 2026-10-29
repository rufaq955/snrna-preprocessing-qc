#!/usr/bin/env python3
"""
Step 01 - Convert 10x Genomics matrix output to AnnData (.h5ad).

Reads the three canonical 10x files (matrix.mtx.gz, features.tsv.gz,
barcodes.tsv.gz) from each sample directory and writes one .h5ad per sample.

Parallelised with a process pool because the gzipped MTX parse is the
bottleneck on multi-sample runs.

Usage
-----
    python scripts/01_10x_to_h5ad.py --samples config/samples.tsv \
                                     --outdir data/interim/h5ad \
                                     --workers 4
"""

from __future__ import annotations

import argparse
import gzip
import logging
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

LOG = logging.getLogger("01_10x_to_h5ad")

REQUIRED_FILES = ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")


def configure_logging(logfile: Path) -> None:
    logfile.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(logfile), logging.StreamHandler(sys.stdout)],
    )


def validate_tsv(path: Path, min_columns: int = 2) -> bool:
    """Cheap sanity check: a TSV row should contain >= min_columns fields.

    Guards against truncated or HTML-error-page downloads, which otherwise
    surface later as a confusing shape mismatch.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", encoding="utf-8") as handle:
            first = handle.readline().strip()
        return len(first.split("\t")) >= min_columns
    except (OSError, UnicodeDecodeError) as exc:
        LOG.warning("Validation failed for %s: %s", path, exc)
        return False


def read_10x_features(path: Path) -> pd.DataFrame:
    """Parse features.tsv.gz keeping BOTH Ensembl ID and gene symbol.

    Column order is (Ensembl ID, gene symbol, feature type). We index by
    symbol because downstream marker inspection and plotting expect
    recognisable gene names, but we retain the Ensembl ID in `.var` so the
    mapping is never lost.
    """
    frame = pd.read_csv(path, sep="\t", header=None, dtype=str)

    if frame.shape[1] >= 2:
        frame = frame.iloc[:, :3] if frame.shape[1] >= 3 else frame
        columns = ["gene_ids", "gene_symbols", "feature_types"][: frame.shape[1]]
        frame.columns = columns
    else:
        frame.columns = ["gene_ids"]

    if "gene_symbols" not in frame.columns:
        frame["gene_symbols"] = frame["gene_ids"]

    # Duplicate symbols are legal in 10x output; make them unique so that
    # AnnData indexing stays well-defined.
    duplicated = frame["gene_symbols"].duplicated(keep=False)
    if duplicated.any():
        LOG.info("Making %d duplicate gene symbols unique", int(duplicated.sum()))
        frame.loc[duplicated, "gene_symbols"] = (
            frame.loc[duplicated, "gene_symbols"]
            + "-"
            + frame.loc[duplicated].groupby("gene_symbols").cumcount().astype(str)
        )

    frame = frame.set_index("gene_symbols", drop=False)
    frame.index.name = None
    return frame


def process_sample(name: str, input_dir: Path, outdir: Path) -> Path | None:
    """Convert a single 10x sample directory into an .h5ad file."""
    paths = {fname: input_dir / fname for fname in REQUIRED_FILES}

    missing = [fname for fname, path in paths.items() if not path.exists()]
    if missing:
        LOG.error("Skipping %s - missing file(s): %s", name, ", ".join(missing))
        return None

    if not validate_tsv(paths["features.tsv.gz"], min_columns=2):
        LOG.error("Skipping %s - features.tsv.gz failed validation", name)
        return None

    try:
        # mtx is genes x cells; transpose to cells x genes and store CSR.
        matrix = scipy.io.mmread(paths["matrix.mtx.gz"]).astype(np.float32)
        matrix = sp.csr_matrix(matrix.T)

        barcodes = pd.read_csv(
            paths["barcodes.tsv.gz"], sep="\t", header=None, dtype=str
        ).iloc[:, 0]
        features = read_10x_features(paths["features.tsv.gz"])

        if matrix.shape != (len(barcodes), features.shape[0]):
            raise ValueError(
                f"Shape mismatch: matrix {matrix.shape} vs "
                f"{len(barcodes)} barcodes x {features.shape[0]} features"
            )

        adata = ad.AnnData(
            X=matrix,
            obs=pd.DataFrame(index=pd.Index(barcodes.astype(str), name="barcode")),
            var=features,
        )
        adata.uns["sample_name"] = name

        outfile = outdir / f"{name}.h5ad"
        adata.write_h5ad(outfile, compression="gzip")
        LOG.info("Wrote %s (%d cells x %d genes)", outfile, *adata.shape)
        return outfile

    except Exception:  # noqa: BLE001 - one bad sample must not kill the run
        LOG.exception("Failed to process sample %s", name)
        return None


def load_samples(path: Path) -> list[tuple[str, Path]]:
    """Read the sample manifest (name<TAB>directory, with a header row)."""
    frame = pd.read_csv(path, sep="\t", comment="#")
    missing_columns = {"name", "directory"} - set(frame.columns)
    if missing_columns:
        raise ValueError(
            f"{path} is missing required column(s): {sorted(missing_columns)}"
        )
    return [(str(r.name), Path(r.directory)) for r in frame.itertuples()]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=Path, default=Path("config/samples.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("data/interim/h5ad"))
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--logfile", type=Path, default=Path("logs/01_convert.log"))
    args = parser.parse_args()

    configure_logging(args.logfile)
    args.outdir.mkdir(parents=True, exist_ok=True)

    samples = load_samples(args.samples)
    LOG.info("Loaded %d sample(s) from %s", len(samples), args.samples)

    failures = 0
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(process_sample, name, directory, args.outdir): name
            for name, directory in samples
        }
        for future in as_completed(futures):
            name = futures[future]
            if future.result() is None:
                failures += 1
                LOG.error("FAILED: %s", name)

    LOG.info("Done. %d succeeded, %d failed.", len(samples) - failures, failures)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
