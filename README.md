# MVAcaller

[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**MVAcaller** is a ready-to-use cell-calling framework for microbial single-cell RNA-seq.
It identifies genuine cells across a wide range of RNA contents by integrating three
complementary statistical tests, and is particularly well suited to datasets where the
classic "knee" of the barcode-rank plot is obscured — the default situation for small
prokaryotic cells, rRNA-dominated libraries and mixed communities such as gut microbiomes.

## Algorithm

MVAcaller (1) fits an ambient-RNA profile from the lowest-UMI barcodes in a user-defined
rank window, (2) computes three P-values — a rRNA-fraction Z-test against the ambient rRNA
percentage, a multinomial EmptyDrops-style test on the KEGG feature space, and a multinomial
test on the non-rRNA gene composition (total RNA in omics mode) — (3) combines them with
Fisher's method, and (4) applies Benjamini-Hochberg correction. Barcodes with adjusted
combined P < `--pvalue-thre` (default 0.01) are called as cells.

## Installation

```bash
git clone https://github.com/MVAcaller.git
cd MVAcaller
pip install -r requirements.txt
```

| Dependency | Needed for |
|---|---|
| `numpy`, `scipy`, `pandas`, `pysam` | cell calling (always) |
| `anndata` | `.h5ad` input / output |
| `bowtie2`, `samtools` | FASTQ → BAM (alignment path only) |
| `seqkit`, `eggNOG-mapper` | building the gene → KEGG map |
| Go 1.19+ | `mva_umi_extract`, `mva_build_sqlite` (optional, faster) |

---

## Two ways to run

Both entry points are shell wrappers that run the whole pipeline in the right order.
They stop with a clear message if a dependency is missing; add `--dry-run` to print every
command without executing it.

```bash
# 1) from paired reads
./run_from_fastq.sh -1 R1.fq.gz -2 R2.fq.gz \
    -n genome.fna -g annotation.gff3 -o out/sample -t 16

# 2) from a count matrix (MTX or .h5ad), e.g. from another aligner
./run_from_matrix.sh --mtx matrix.mtx.gz --barcodes barcodes.tsv.gz \
    --features features.tsv.gz \
    -n genome.fna -g annotation.gff3 -o out/sample -t 16
```

`genome.fna` and `annotation.gff3` are the files downloaded from the assembly database
(e.g. NCBI Datasets). The gene → KEGG map is built from them automatically
(`mva_extract_features.py` → `seqkit translate` → `emapper.py`), so no extra input is
required. Options that affect that build:

| Flag | Meaning |
|---|---|
| `--kegg-map FILE` | reuse a map built earlier instead of rebuilding it |
| `--eggnog-db DIR` | eggNOG data directory; defaults to the one your eggNOG-mapper is configured with |

The map is cached in `out/sample/ref/gene_module.tsv` (FASTQ path) or
`out/sample/work/gene_module.tsv` (matrix path) and reused on later runs.

### Order and logic

**From FASTQ (`run_from_fastq.sh`)**

1. Annotation → GFF3 (a GTF is converted first). The GFF3 is only used to recognise rRNA genes.
2. Genome + GFF3 → gene-tagged FASTA (`cds` / `trna` / `rrna` / `ncrna`) and a Bowtie2 index.
3. Gene → KEGG map: translate the CDS from step 2 and annotate them with eggNOG-mapper.
   The resulting `gene_module.tsv` is cached and reused; `--kegg-map` skips this step.
4. Move barcode and UMI from R1 into the R2 read name.
5. Align R2 to the gene-tagged index, so each alignment maps to a gene rather than a genome.
6. Deduplicate UMI: one molecule per (barcode, gene, UMI).
7. Build a SQLite molecule database, then run MVAcaller on it.
8. Write the filtered gene × cell matrix for the barcodes called as cells.

**From a matrix (`run_from_matrix.sh`)**

1. Annotation → GFF3 (same as above).
2. Gene → KEGG map, built from `--genome` + the GFF3 unless `--kegg-map` is given.
3. MTX / `.h5ad` → the same SQLite molecule database.
4. MVAcaller on that database.

No alignment and no Bowtie2 index are involved in the matrix path. The matrix **must still
contain the ambient (empty-droplet) barcodes**: the ambient profile is learned from the
low-UMI tail of the barcode rank curve, and a cell-filtered matrix silently produces a
biased profile.

### Platform formats are not uniform

Every microbial scRNA-seq platform places the barcode and UMI differently, so step 4 is the
only step that cannot be generic. `run_from_fastq.sh` exposes it as `--bc-start`, `--bc-len`,
`--umi-start`, `--umi-len`. If your platform does not fit (pre-extracted barcodes, custom
read names, CB/UB tags), pass `--skip-extract`, bring your own gene-tagged BAM, and set the
matching MVAcaller options (`--qname-regex`, `--use-cb-tag`, `--gene-delimiter`).

### Running MVAcaller alone

If you already have a deduplicated gene-tagged BAM or a SQLite database:

```bash
# first run builds molecules.sqlite, later runs reuse it
python MVAcaller.py -i dedup.bam -g annotation.gff3 --kegg-map gene_module.tsv \
    --build-sqlite molecules.sqlite --reuse-sqlite -o out/mva_run --threads 16

# reuse only (no -i needed)
python MVAcaller.py -g annotation.gff3 --kegg-map gene_module.tsv \
    --build-sqlite molecules.sqlite --reuse-sqlite -o out/mva_run --threads 16
```

Add `--mode-omics` for multi-species datasets (human / mouse gut), which replaces the
non-rRNA composition test with a total-RNA test.

## Reference files

The three reference inputs are built once per genome / species set. Both wrapper scripts
do this for you; the steps below are what they run, in order, if you prefer to do it by hand:

1. `mva_extract_features.py` — splits the annotation into CDS / tRNA / rRNA / ncRNA and
   writes `features.fa` with `><seqid>#<gene_id>` headers, the RNAME convention MVAcaller
   expects (mobile elements are skipped unless `--keep-mobile`).
2. `seqkit translate ref/cds.fa > protein.faa` — protein sequences for annotation.
3. `emapper.py -i protein.faa --itype proteins --data_dir <eggnog_db> -o eggnog_out` —
   functional annotation.
4. `grep -v '^#' eggnog_out.emapper.annotations | cut -f1,14 | sed 's/^[^#]*#//' > gene_module.tsv`
   — column 14 is the KEGG module; the `sed` removes the `<contig>#` prefix so IDs match
   the bare gene IDs used by MVAcaller.
5. `bowtie2-build ref/features.fa ref/features` — the gene-tagged alignment index.
6. `mva_dedup.py` — UMI deduplication; it also writes `counts.tsv` (barcode, reads, UMIs)
   for `mva_knee.py` / `mva_rank_target.py`.

For multi-species work, concatenate the per-species GFF3s and build one `features.fa` and
one index for the whole community.

`--gff` must be GFF3 and rRNA genes must carry `gene_biotype=rRNA` (Prokka/Bakta style).
If you only have a GTF, convert it with `mva_gtf2gff3.py`.

## Input formats

| Input | Supported | Note |
|---|---|---|
| Gene-tagged BAM / SAM (`.sam.gz`) | native | in-memory or via SQLite |
| SQLite molecule DB | native | `--build-sqlite ... --reuse-sqlite`, no `-i` |
| Raw MTX triplets, `.h5ad`, molecule table | via `mva_mtx2sqlite.py` | matrix must include empty droplets |
| Cell-filtered matrix | runs but biased | no ambient barcodes left |
| `.loom`, Seurat | no | export to MTX / `.h5ad` first |

## Output

| File | Content |
|---|---|
| `out/mva_run.cells.txt` | barcodes called as cells, one per line |
| `out/mva_run.barcodes.tsv` | per-barcode table: UMI counts, rRNA fraction, three P-values, adjusted and combined P, `call`, `reason` |
| `out/mva_run.summary.txt` | read / cell counts, ambient window and thresholds actually used |
| `out/mva_run.filtered_matrix.sparse.tsv` | gene × cell counts restricted to the called cells, sparse COO (`gene`, `cell`, `count`); add `--dense-matrix` for a dense table instead |

`barcode`, `call` and `reason` are the columns needed to hand the result to Scanpy or
Seurat; `reason` records why a barcode was retained or rejected (`combined_padj`,
`umi_retained`, `above_upper_bound`, ...).

The sparse matrix is written from the deduplicated BAM, so each entry is one molecule.
To continue in Scanpy, convert it to AnnData:

```bash
python mva_load_h5ad.py -i out/mva_run.filtered_matrix.sparse.tsv -o out/mva_run.filtered.h5ad
```

## Key parameters

| Flag | Default | Meaning |
|---|---|---|
| `--pvalue-thre` | 0.01 | combined P threshold after BH correction |
| `--ind-min` / `--ind-max` | auto Q3 rank / 570000000 | ambient window (rank range) |
| `--umi-lower` | 300 (100 in omics) | barcodes at or below are skipped as empty |
| `--umi-upper` | 10000 | barcodes at or above are called without a test |
| `--soft-umi-thre` | on | UMI-dependent P threshold; off via `--no-soft-umi-thre` |
| `--sim-n` | 10000 | multinomial simulations per barcode |
| `--mode-omics` | off | total-RNA third test for multi-species data |
| `--dense-matrix` | off | wrapper flag: dense gene × cell table instead of sparse COO |
| `--threads` | 1 | workers for the simulations |

`python MVAcaller.py --help` lists every option.

## Repo layout

```
MVAcaller/
├── MVAcaller.py              main cell-calling program
├── run_from_fastq.sh         pipeline: FASTQ + genome + annotation → cells
├── run_from_matrix.sh        pipeline: MTX / .h5ad + annotation → cells
├── mva_umi_extract.go        barcode / UMI from R1 into the R2 read name
├── mva_extract_features.py   CDS / tRNA / rRNA / ncRNA + gene-tagged FASTA
├── mva_gtf2gff3.py           GTF → GFF3
├── mva_dedup.py              UMI deduplication → dedup.bam, counts.tsv
├── mva_build_sqlite.go       BAM / SAM → SQLite molecule DB
├── mva_mtx2sqlite.py         MTX / .h5ad / molecule table → SQLite
├── mva_knee.py               knee / inflection point on counts.tsv
├── mva_rank_target.py        rank-based barcode subsetting
├── mva_load_h5ad.py          sparse TSV → AnnData for Scanpy
├── mva_bam_to_mat.py         dedup BAM + called barcodes → filtered gene × cell matrix
├── benchmarks/               paper-reproduction scripts
├── requirements.txt
└── LICENSE
```

## Citation


## License

MIT — see [LICENSE](LICENSE).
