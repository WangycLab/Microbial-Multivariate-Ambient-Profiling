# MVAcaller

[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**MVAcaller** is a ready-to-use cell-calling framework for microbial single-cell RNA-seq.
It identifies genuine cells across a wide range of RNA contents by integrating three
complementary statistical tests and is particularly well suited to datasets where the
classic "knee" of the barcode-rank plot is obscured — which is the default situation for
small prokaryotic cells, rRNA-dominated libraries, and complex mixed communities such as
gut microbiomes.

## Algorithm

For each barcode, MVAcaller (1) fits an ambient-RNA profile from the lowest-UMI barcodes
user-defined rank window, (2) computes three P-values: a rRNA-fraction Z-test against the
ambient rRNA percentage, a multinomial EmptyDrops-style test on the KEGG-KO feature space,
and a multinomial test on the non-rRNA gene composition (or total-RNA in multi-species /
omics mode), (3) combines the three P-values with Fisher's method, and (4) applies
Benjamini-Hochberg correction. Barcodes with adjusted combined P < `--pvalue-thre`
(default **0.01**) are called as cells.

## Installation

```bash
git clone https://github.com/MVAcaller.git
cd MVAcaller
pip install -r requirements.txt
```
### Dependencies

| Package | Purpose |
|---|---|
| `numpy`, `scipy`, `pandas` | Core stats |
| `pysam` | BAM/SAM parsing (required if your input is BAM/CRAM) |
| `anndata`, `scikit-learn` | Optional, only for `mmd_load_h5ad.py` |
| `bowtie2`, `samtools` | Optional, for the preprocessing utilities only |

---



