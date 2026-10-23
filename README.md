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
| `seqkit`, `eggNOG-mapper` | Optional, only for the one-time reference build |

---

## Building the reference (one-time prep)

MVAcaller needs three reference files that are built once per genome / species set:
a **gene-tagged alignment index**, a **GFF3 annotation** (for rRNA genes), and a
**gene → KEGG module (KM) map**.

```bash
python mva_extract_features.py -g annotation.gtf -f genome.fna -o ref/
```

Outputs `ref/cds.fa`, `ref/trna.fa`, `ref/rrna.fa`, `ref/ncrna.fa` and the combined
`ref/features.fa`. Every header is written as `><seqid>#<gene_id>`, which is exactly the
RNAME convention MVAcaller expects. Mobile elements (insertion sequences, prophages,
plasmids, integrons) are skipped by default; add `--keep-mobile` to retain them.

```bash
seqkit translate ref/cds.fa > protein.faa
```

Use `seqkit translate -f 11` (bacterial/archaeal code) if your organisms use
alternative start codons beyond the standard table.

### Step 3 — functional annotation with eggNOG-mapper

```bash
emapper.py -i protein.faa --itype proteins \
    --data_dir /path/to/eggnog_db \
    -o eggnog_out --cpu 64 --dbmem
```

`--data_dir` points at a downloaded eggNOG database (`emapper.py --download_eggnog_data.py`).
`--dbmem` loads the database into RAM (faster, needs ~30–60 GB for the full set).
This produces `eggnog_out.emapper.annotations`.

### Step 4 — build `gene_module.tsv`

```bash
grep -v '^#' eggnog_out.emapper.annotations \
  | cut -f1,14 \
  | sed 's/^[^#]*#//' \
  > gene_module.tsv
```
### Step 5 — build the gene-tagged alignment index

```bash
bowtie2-build ref/features.fa ref/features
```

Reads must then be aligned to this index, so that each RNAME is `<contig>#<gene_id>`:

```bash
bowtie2 -x ref/features -U reads.fq --very-sensitive-local -p 16 \
  | samtools view -bS - \
  | samtools sort -@ 8 -o temp.sorted.bam
```

### Step 6 — UMI deduplication

```bash
python mva_dedup.py -b temp.sorted.bam -c counts.tsv -d dedup.bam -o dedup.fq
```

One pass over the sorted BAM keeps, per barcode, the first read seen for every
**(contig, gene, UMI)** combination

| Output | Contents |
|---|---|
| `dedup.bam` | deduplicated gene-tagged BAM — **this is the input for `MVAcaller.py`** |
| `dedup.fq`  | the retained reads as FASTQ (for re-alignment, ambient simulation, etc.) |
| `counts.tsv` | header + `barcode`, `reads_count`, `umi_count`, sorted by UMI count descending |

### Step 7 — annotation for `--gff`

`--gff` must be **GFF3** (it is used to pick out rRNA genes). If you only have a GTF:

```bash
python mva_gtf2gff3.py annotation.gtf annotation.gff3
```

For multi-species datasets, concatenate (cat) the per-species GFF3s and build a single
`features.fa` / gene-tagged index for the whole community.

---

### Required input formats

`MVAcaller.py` performs cell calling **directly from aligned reads**, so its input needs to be
a **gene-tagged BAM / SAM** (optionally `.sam.gz`).

# We recommend to first build SQlite indix for running MVA, which can be reused if necessarily, instead of uploaded to memory once and released after running, via

go build -o mva_build_sqlite mva_build_sqlite.go
./mva_build_sqlite -i dedup.bam -o molecules.sqlite

# We realised users may obtain the gene expression matrix or h5ad through different ways, thus an alternative approach to build index from mtx/h5ad to SQlite is provided, via:

# 10x-style triplets (raw_feature_bc_matrix)
python mva_mtx2sqlite.py --mtx matrix.mtx.gz --barcodes barcodes.tsv.gz \
    --features features.tsv.gz --strip-suffix=-1 -o molecules.sqlite
    
# AnnData from another aligner
python mva_mtx2sqlite.py --h5ad raw.h5ad --strip-suffix=-1 -o molecules.sqlite

# cell calling on the resulting database
python MVAcaller.py --build-sqlite sample.sqlite --reuse-sqlite \
    -i sample.sqlite -g annotation.gff3 --kegg-map gene_module.tsv -o out/sample


| Flag | Description | Format |
|---|---|---|
| `-i, --input` | Gene-tagged BAM (sorted) or SAM (optionally `.sam.gz`) | BAM / SAM / SAM.gz |
| `-g, --gff` | GFF3 annotation used to identify rRNA genes | GFF3 |
| `--kegg-map` | **Two columns:** `gene_id <tab> KO_number` | plain TSV, no header |

### Optional inputs

| Flag | Description |
|---|---|
| `--build-sqlite PATH` | Build and (on rerun) reuse a SQLite molecule DB. Recommended for BAMs > ~500 MB. |
| `--target-barcodes FILE` | Restrict testing to barcodes in this file (one per line) — useful for fast iteration on a subset. |






