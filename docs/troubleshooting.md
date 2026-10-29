# Troubleshooting

Notes on the failure modes that actually come up with this stack.

## Step 01 — conversion

**`Shape mismatch: matrix (N, M) vs X barcodes x Y features`**
The three 10x files are out of sync — usually a truncated download or a
directory holding output from a different run. Re-download from the source
rather than trying to patch it.

**Gene symbols are missing / `var_names` are Ensembl IDs**
You are reading `features.tsv.gz` with only column one. Step 01 keeps both
`gene_ids` and `gene_symbols`; check `adata.var.head()`.

**Memory blow-up on large samples**
`scipy.io.mmread` materialises the dense intermediate. For very large samples,
reduce `--workers` so fewer matrices are held in RAM simultaneously.

## Step 02 — CellBender

**`cellbender: command not found`**
CellBender is usually installed inside a conda env. Either activate the env or
pass an explicit path: `--cellbender /path/to/envs/<env>/bin/cellbender`.

**CUDA errors**
Drop `--cuda` and retry on CPU. If the CUDA build of PyTorch and the driver
disagree, CellBender fails at model init rather than gracefully falling back.

**Output is degenerate (near-zero counts)**
Usually too few cells or too little depth for the model to separate cells from
ambient. CellBender is not meaningful on toy datasets — test steps 01/03/06
with synthetic data and 02/04/05 on real samples.

## Step 03 — ptrepack

**`ptrepack not found`**
`pip install tables`.

**`Seurat cannot open the h5 file`**
The `:/matrix` group suffix matters. Copying the whole file rather than the
`/matrix` group produces a layout `Read10X_h5` will not accept.

## Step 04 — DoubletFinder

**`paramSweep` is very slow**
Cost scales with cell number and PC range. For samples above ~20k cells,
consider reducing `--pcs`.

**Doublet counts look implausibly high**
`--doublet-rate` is an assumption, not a measurement. 7.5% suits roughly
5,000 recovered cells; higher loading rates need a higher value. Cross-check
against expected multiplet rate for your 10x chip.

## Step 05 — QC

**Everything is filtered out**
`percent.mt` computed as `NA` or `0` across the board usually means the gene
symbols are Ensembl IDs, so the `^MT-` pattern matches nothing. Confirm
`var_names` are symbols (see step 01).

**Nuclei pass with suspiciously high mito fraction**
This is snRNA-seq, not scRNA-seq — nuclei should show low mitochondrial
content. If the distribution looks like whole cells, the input may be cells
rather than nuclei; adjust `--max-mt` accordingly.

## Step 06 — h5ad ↔ rds

**`SeuratDisk::Convert` fails on a v5 assay**
SeuratDisk lags behind Seurat's v5 assay changes. Load with
`LoadH5Seurat(..., assays = "RNA")`, or downgrade the object to a v3 assay
with `CreateSeuratObject` before converting.

**`hdf5r` fails to install**
Needs HDF5 headers. On conda: `conda install -c conda-forge hdf5 r-hdf5r`.

**Metadata attached but cell names silently reordered**
Step 06b asserts that cell names match, and reorders only when the sets are
equal. If it warns instead of attaching, the metadata CSV came from a
different object — re-run step 06a against the same `.h5ad`.
