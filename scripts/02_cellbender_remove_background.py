#!/usr/bin/env python3
"""
Step 02 - Ambient RNA removal with CellBender.

Runs `cellbender remove-background` over every .h5ad produced by step 01.
CellBender emits `*_cellbender_filtered.h5`, which step 03 repacks for R.

Samples whose output already exists are skipped, so the step is safe to
re-run after a partial or interrupted run.

Usage
-----
    python scripts/02_cellbender_remove_background.py \
        --indir data/interim/h5ad \
        --outdir data/interim/cellbender \
        --cuda
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

LOG = logging.getLogger("02_cellbender")

# Defaults mirror the parameters used for the snRNA-seq cohort. They are
# exposed on the command line so a run can be reproduced or tuned without
# editing this file.
DEFAULT_CELLBENDER_ARGS = [
    "--projected-ambient-count-threshold", "5",
    "--checkpoint-mins", "120",
    "--learning-rate", "0.0001",
    "--training-fraction", "0.9",
    "--low-count-threshold", "20",
    "--epochs", "150",
    "--posterior-batch-size", "200",
]


def configure_logging(logfile: Path) -> None:
    logfile.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(logfile), logging.StreamHandler(sys.stdout)],
    )


def resolve_cellbender(explicit: str | None) -> str:
    """Find the cellbender executable.

    CellBender is often installed inside a conda env, so we accept an
    explicit path before falling back to PATH lookup.
    """
    if explicit:
        path = Path(explicit)
        if not path.exists():
            raise FileNotFoundError(f"--cellbender path does not exist: {explicit}")
        return str(path)

    found = shutil.which("cellbender")
    if not found:
        raise FileNotFoundError(
            "cellbender not found on PATH. Install it (see environment.yml) "
            "or pass --cellbender /path/to/cellbender."
        )
    return found


def run_cellbender(
    cellbender: str,
    input_file: Path,
    output_prefix: Path,
    extra_args: list[str],
) -> tuple[bool, str]:
    """Execute remove-background for one sample. Returns (ok, stderr_tail)."""
    command = [
        cellbender, "remove-background",
        "--input", str(input_file),
        "--output", str(output_prefix),
        *extra_args,
    ]

    LOG.info("Running: %s", " ".join(command))
    started = time.time()

    completed = subprocess.run(
        command, capture_output=True, text=True, check=False
    )
    elapsed = time.time() - started

    if completed.stdout:
        LOG.debug("STDOUT:\n%s", completed.stdout)
    if completed.stderr:
        LOG.debug("STDERR:\n%s", completed.stderr)

    if completed.returncode == 0:
        LOG.info("Completed %s in %.1fs", input_file.name, elapsed)
        return True, ""

    LOG.error(
        "CellBender failed for %s (exit %d)", input_file.name, completed.returncode
    )
    return False, completed.stderr[-2000:]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indir", type=Path, default=Path("data/interim/h5ad"))
    parser.add_argument("--outdir", type=Path, default=Path("data/interim/cellbender"))
    parser.add_argument("--logfile", type=Path, default=Path("logs/02_cellbender.log"))
    parser.add_argument("--cellbender", default=None, help="Explicit path to cellbender")
    parser.add_argument("--cuda", action="store_true", help="Add --cuda")
    parser.add_argument(
        "--cellbender-arg", action="append", default=[],
        metavar="ARG",
        help="Override/extend CellBender args (repeatable), e.g. --cellbender-arg=--epochs --cellbender-arg=100",
    )
    args = parser.parse_args()

    configure_logging(args.logfile)

    if not args.indir.is_dir():
        LOG.error("Input directory does not exist: %s", args.indir)
        return 2

    cellbender = resolve_cellbender(args.cellbender)
    args.outdir.mkdir(parents=True, exist_ok=True)

    extra_args = args.cellbender_arg or DEFAULT_CELLBENDER_ARGS
    if args.cuda:
        extra_args = [*extra_args, "--cuda"]

    h5ad_files = sorted(args.indir.glob("*.h5ad"))
    if not h5ad_files:
        LOG.error("No .h5ad files found in %s", args.indir)
        return 2

    LOG.info("Found %d sample(s); started %s", len(h5ad_files), datetime.now())

    succeeded, skipped, failed = 0, 0, []

    for index, input_file in enumerate(h5ad_files, start=1):
        prefix = args.outdir / f"{input_file.stem}_cellbender"
        expected = args.outdir / f"{input_file.stem}_cellbender_filtered.h5"

        if expected.exists():
            LOG.info("[%d/%d] Already done, skipping: %s", index, len(h5ad_files), input_file.name)
            skipped += 1
            continue

        LOG.info("[%d/%d] Processing: %s", index, len(h5ad_files), input_file.name)
        try:
            ok, error_tail = run_cellbender(cellbender, input_file, prefix, extra_args)
        except KeyboardInterrupt:
            LOG.warning("Interrupted by user; stopping.")
            break
        except Exception as exc:  # noqa: BLE001
            ok, error_tail = False, str(exc)

        if ok:
            succeeded += 1
        else:
            failed.append((input_file.name, error_tail))

    LOG.info("=" * 60)
    LOG.info("Total: %d | succeeded: %d | skipped: %d | failed: %d",
             len(h5ad_files), succeeded, skipped, len(failed))
    for name, error in failed:
        LOG.error("  FAILED %s\n%s", name, error)

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
