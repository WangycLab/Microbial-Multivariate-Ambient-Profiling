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
| `numpy`, `scipy`, `pandas` | Core statistics |
| `pysam` | BAM/SAM parsing |
| `anndata` | Optional, only for `.h5ad` input/output (`mva_mtx2sqlite.py`, `mva_load_h5ad.py`) |
| `samtools` | Optional, but required when the Go builder reads a BAM |
| `bowtie2` | Optional, only for the preprocessing utilities |
| `seqkit`, `eggNOG-mapper` | Optional, only for the one-time reference build |
| Go 1.19+ | Optional, to build `mva_build_sqlite` and `mva_umi_extract` |

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

### Step 1 — functional annotation with eggNOG-mapper

```bash
emapper.py -i protein.faa --itype proteins \
    --data_dir /path/to/eggnog_db \
    -o eggnog_out --cpu 64 --dbmem
```

`--data_dir` points at a downloaded eggNOG database (`emapper.py --download_eggnog_data.py`).
`--dbmem` loads the database into RAM (faster, needs ~30–60 GB for the full set).
This produces `eggnog_out.emapper.annotations`.

### Step 2 — build `gene_module.tsv` (KM, kegg module)

```bash
grep -v '^#' eggnog_out.emapper.annotations \
  | cut -f1,14 \
  | sed 's/^[^#]*#//' \
  > gene_module.tsv
```
### Step 3 — build the gene-tagged alignment index

```bash
bowtie2-build ref/features.fa ref/features
```

Reads must then be aligned to this index, so that each RNAME is `<contig>#<gene_id>`:

```bash
bowtie2 -x ref/features -U reads.fq --very-sensitive-local -p 16 \
  | samtools view -bS - \
  | samtools sort -@ 8 -o temp.sorted.bam
```

### Step 4 — UMI deduplication

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

### Step 5 — annotation for `--gff`

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

**Build the SQLite DB once and reuse it.** It is much faster than reloading a BAM for every
threshold change, and it is the only input needed on reruns:

```bash
go build -o mva_build_sqlite mva_build_sqlite.go  
./mva_build_sqlite -i dedup.bam -o molecules.sqlite 
```

Users who obtain a gene expression matrix from another aligner can build the identical
database from MTX / `.h5ad` instead:

```bash
# 10x-style triplets (raw_feature_bc_matrix)
python mva_matrix_to_sqlite.py --mtx matrix.mtx.gz --barcodes barcodes.tsv.gz \
    --features features.tsv.gz --strip-suffix=-1 -o molecules.sqlite

# AnnData from another aligner
python mva_matrix_to_sqlite.py --h5ad raw.h5ad --strip-suffix=-1 -o molecules.sqlite
```

## Quick start

### 1. from SQlite built (recommended)

```bash
python MVAcaller.py \
    -g annotation.gff3 \
    --kegg-map gene_module.tsv \
    --build-sqlite molecules.sqlite \
    --reuse-sqlite \
    -o out/mva_run \
    --threads 16
```

### 2. In-memory mode, from bam

```bash
python MVAcaller.py \
    -i dedup.bam \
    -g annotation.gff3 \
    --kegg-map gene_module.tsv \
    -o out/mva_run \
    --threads 16
```

### 3. Multi-species / omics mode (e.g., multi-species/ human/ mouse gut microbiomes)

```bash
python MVAcaller.py \
    --build-sqlite molecules.sqlite \
    --reuse-sqlite \
    -g annotations_combined.gff3 \
    --kegg-map gene_module.tsv \
    -o out/mva_run \
    --mode-omics \
    --threads 24
```

## Output files

For every run with `-o out/mva_run`, MVAcaller writes:

| File | Content |
|---|---|
| `out/mva_run.barcodes.tsv` | One row per tested barcode (full P-value table, see below) |
| `out/mva_run.cells.txt` | Just the barcode IDs called as cells — one per line |
| `out/mva_run.summary.txt` | Read counts, cell counts, ambient-profile settings, thresholds used |

### Columns in `*.barcodes.tsv`

| Column | Description |
|---|---|
| `barcode` | Cell barcode |
| `rank` | Rank by total UMI (1 = highest) |
| `umis` | Total UMI count |
| `rrna_umis` / `rrna_pct` | UMI count and fraction assigned to rRNA genes |
| `ko_umis` / `non_rrna_umis` / `non_rrna_genes` | KEGG-KO, non-rRNA UMI counts and non-rRNA gene count |
| `genes` | Total number of genes detected |
| `ko_features` | Number of KEGG-KM features detected |
| `call` | **`cell`** or **`empty`** |
| `reason` | Why: `combined_padj`, `umi_retained`, `above_upper_bound`, `below_min_non_rrna_genes`, etc. |
| `pvalue_rrnapct`, `pvalue_ko`, `pvalue_non_rrna` | Three per-test raw P-values |
| `padj_rrnapct`, `padj_ko`, `padj_non_rrna` | Three per-test BH-adjusted P-values |
| `pvalue_combined`, `padj_combined` | Fisher-combined raw and BH-adjusted P-values |
| `pvalue_thre_used` | P-threshold actually applied (equals `--pvalue-thre`, or the interpolated one under `--soft-umi-thre`) |

## Parameters

### Ambient profile

| Flag | Default | Description |
|---|---|---|
| `--ind-min` | auto | Start rank (inclusive) of the ambient window. Auto = Q3 rank, i.e. the bottom 25% of barcodes by UMI |
| `--ind-max` | 570000000 | End rank (inclusive); effectively the end of the tail |
| `--ambient-pseudocount` | 1.0 | Pseudocount added to the ambient profile |

### Statistics

| Flag | Default | Description |
|---|---|---|
| `--pvalue-thre` | 0.01 | Combined P threshold (post-BH) for calling cells |
| `--sim-n` | 10000 | Multinomial simulations per barcode |
| `--threads` | 1 | Worker threads for the simulations |
| `--seed` | 19760110 | Random seed |
| `--no-bh` | off | Call cells on the raw Fisher P instead of the BH-adjusted value |
| `--mode-omics` | off | Use a total-RNA test as the third test; `--umi-lower` defaults to 100 |

### UMI thresholds

| Flag | Default | Description |
|---|---|---|
| `--umi-lower` | 300 (100 with `--mode-omics`) | Barcodes with UMI ≤ this are skipped as empty |
| `--umi-upper` | 10000 | Barcodes with UMI ≥ this are called as cells without computing P |
| `--umi-thre` | off | Barcodes above this are forced to cell (P = 0) but still enter BH correction |
| `--min-non-rrna-genes` | 0 (off) | Skip barcodes with fewer detected non-rRNA genes |
| `--soft-umi-thre` | **on** | Slide the P threshold along log(UMI): stricter for low-UMI barcodes. Disable with `--no-soft-umi-thre` |
| `--soft-umi-low` / `--soft-umi-high` | 100 / 1000 | UMI anchors between which the threshold is interpolated |
| `--soft-thre-low` / `--soft-thre-high` | 1e-5 / 0.01 | Threshold at the low and high anchors |

### BAM/SAM parsing

| Flag | Default | Description |
|---|---|---|
| `--qname-simple` | **on** | Parse `..._<BC>_<UMI>` from the read name. Disable with `--no-qname-simple` |
| `--qname-regex` | `_(?P<cb>[ACGTN]+)_(?P<umi>[ACGTN]+)$` | Regex used when `--qname-simple` is off |
| `--use-cb-tag` | off | Read CB/UMI from SAM tags (`CB:Z:` plus `UB:Z:` / `UR:Z:`) instead of the name |
| `--gene-delimiter` | `#` | Delimiter between contig and gene ID in RNAME |
| `--min-mapq` | 0 | Minimum MAPQ |
| `--keep-secondary` | off | Keep secondary / supplementary alignments |

### Other

| Flag | Default | Description |
|---|---|---|
| `--multi-ko` | `split` | Genes annotated to several KOs: `split` (counts split equally), `first`, `drop` |
| `--batch-size` | 200000 | SQLite insert batch size |
| `--progress-every` | 1000000 | Progress line every N reads; 0 disables |
| `--log-likelihood` | off | Write per-barcode log-multinomial likelihood to `<prefix>.loglik.tsv` |

## Repo layout

```
MVAcaller/
├── MVAcaller.py                main cell-calling program
├── mva_mtx2sqlite.py     raw MTX / .h5ad / molecule table → SQLite
├── mva_build_sqlite.go         fast BAM/SAM → SQLite builder (optional, Go)
├── mva_umi_extract.go          R1 BC/UMI → R2 header extraction (optional, Go)
├── mva_extract_features.py     CDS / tRNA / rRNA / ncRNA extraction + gene-tagged FASTA
├── mva_dedup.py                UMI deduplication (BAM → dedup.bam / .fq / counts.tsv)
├── mva_gtf2gff3.py          GTF → GFF3 conversion
├── mva_knee.py                 knee / inflection point on counts.tsv
├── mva_rank_target.py          rank-based barcode subsetting and UMI targets
├── mva_load_h5ad.py            sparse TSV → AnnData (.h5ad) for Scanpy
├── mva_bam_to_mat.py           gene-tagged BAM → sparse count matrix
├── benchmarks/                 paper-reproduction scripts on simulated datasets
├── requirements.txt
└── LICENSE
```




