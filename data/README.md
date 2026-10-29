# Data

**This repository contains no data.** It is a code-only pipeline.

Raw sequencing output and processed count matrices from human samples cannot
be redistributed, and de-identified matrices are often still covered by
institutional data-use agreements. Nothing in `data/` is tracked by git — see
`.gitignore`.

## Running on your own data

Point `config/samples.tsv` at directories containing the three canonical 10x
files:

```
your_sample/
├── matrix.mtx.gz
├── features.tsv.gz
└── barcodes.tsv.gz
```

```bash
cp config/samples.example.tsv config/samples.tsv
# edit the paths, then:
bash run_pipeline.sh --samples config/samples.tsv --outdir data
```

## Generating a synthetic dataset for testing

To verify the pipeline runs end-to-end without real data, create a small
random 10x-format sample:

```bash
python - <<'PY'
import gzip, numpy as np, scipy.io, scipy.sparse
from pathlib import Path

out = Path("data/synthetic/sample_01"); out.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(0)

n_genes, n_cells = 2000, 300
density = 0.08
mask = rng.random((n_genes, n_cells)) < density
counts = scipy.sparse.coo_matrix(
    (rng.poisson(2, mask.sum()), np.where(mask)),
    shape=(n_genes, n_cells),
)
scipy.io.mmwrite(out / "matrix.mtx", counts)
with gzip.open(out / "matrix.mtx.gz", "wb") as f:
    f.write((out / "matrix.mtx").read_bytes())
(out / "matrix.mtx").unlink()

genes = [f"GENE{i:05d}" for i in range(n_genes)]
with gzip.open(out / "features.tsv.gz", "wt") as f:
    for i, g in enumerate(genes):
        f.write(f"ENSG{i:011d}\t{g}\tGene Expression\n")

with gzip.open(out / "barcodes.tsv.gz", "wt") as f:
    for i in range(n_cells):
        f.write(f"CELL{i:05d}-1\n")

print("Wrote synthetic sample to", out)
PY
```

Then add a matching row to `config/samples.tsv`:

```
name	directory
synthetic_01	data/synthetic/sample_01
```

Note that step 02 (CellBender) needs enough cells and depth to fit its model;
a 300-cell toy sample may fail or return degenerate output. Use the synthetic
data to validate steps 01, 03 and 06, and test 02/04/05 on a real sample.
