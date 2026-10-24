# benchmarks/ — synthetic-mixture benchmarking

This directory reproduces the benchmarking in the MVAcaller manuscript: each replicate
builds a mixture of **large cells**, **downsampled ("small") cells** and **simulated empty
droplets**, aligns it to the gene-tagged reference, and compares MVAcaller's calls against a
knee-based UMI threshold. Because every barcode's class is known by construction, recall and
FDR can be computed exactly.

These scripts are not part of the end-user toolchain. For calling cells on real data use
[`run_from_fastq.sh`](../run_from_fastq.sh) or [`run_from_matrix.sh`](../run_from_matrix.sh)
— see the main [README](../README.md).

## Directory contents

| File | Purpose |
|---|---|
| `run_benchmark.sh` | Per-sample wrapper: reference build once, then `N` replicates of mix → align → call → summarise. |
| `summarize_calls.awk` | Counts called barcodes by class from an MVAcaller `*.barcodes.tsv`. |
| `tools/split_by_rank_index.py` | Splits a deduplicated FASTQ by barcode UMI rank into large-cell and ambient reads. |
| `tools/simulate_empty2.py` | Resamples reads from the ambient FASTQ into artificial `EMPTY*` barcodes. |
| `tools/gene_proportional_sampley.py` | Downsamples real cell reads in a gene/KEGG-module-proportional way; barcodes get a `#small` suffix. |
| `tools/count_bc.py` | `barcode <tab> reads <tab> UMI` table of a FASTQ (for the knee method). |

Repository-root helper used here: [`mva_knee.py`](../mva_knee.py) (kneedle threshold).

## Benchmark construction

The **large-cell and ambient read sets are built once per sample**; only the empty and
small-cell mixtures are re-drawn between replicates.

Reference build (once):

1. `tools/split_by_rank_index.py` ranks barcodes by read count and splits `dedup.fq`
   into the top `--large-n` barcodes (`<sample>.large.fq`) and everything ranked after
   `--min-ind` (`<sample>.ambient.fq`); reads in between are dropped;
2. `sed 's/\(_[ACGT]*\)_/\1#large_/'` appends `#large` to the large barcodes
   (`<sample>.large.labeled.fq`);
3. the large cells are aligned (`bowtie2 --very-sensitive-local`) and sorted
   (`<sample>.large.sorted.bam`); the same index is used for every mixture.

Then, for each replicate (default 3):

4. `tools/simulate_empty2.py` draws `--empty-n` empty droplets from `<sample>.ambient.fq`,
   mean `--empty-mean` reads per barcode, barcodes prefixed `EMPTY` (seeds 51/52/53);
5. `tools/gene_proportional_sampley.py` downsamples the large cells with coefficients in
   [`--coef-min`, `--coef-max`], producing `--small-n` barcodes suffixed `#small`
   (seeds 51/52/53);
6. empty + small + large reads are concatenated, aligned to the same index and sorted;
7. MVAcaller runs with `--mode-omics`, `--umi-thre 10000`, `--sim-n 10000`,
   `--ind-min 300000 --ind-max 865739999` (seeds 57/58/59);
8. per class, barcodes called by the combined test and by the non-rRNA test are counted,
   and the knee UMI threshold (`mva_knee.py --sensitivity 0.7`) is applied to the same
   mixed library for comparison.

Barcode class is encoded in the name, so the labels are unambiguous:
`EMPTY*` = empty droplet, `*#small` = downsampled cell, `*#large` and any unmarked
barcode = large cell.

## Inputs you must provide

| Flag | Req. | What it is |
|---|---|---|
| `--dedup-fq` | yes | Deduplicated FASTQ with `..._<BC>_<UMI>` in the read name (`dedup.fq` from `mva_dedup.py`) |
| `--gff` | yes | GFF3 annotation (`mva_gtf2gff3.py`, or native GFF3) |
| `--kegg-map` | yes | `gene <tab> KEGG_module`; build it with eggNOG-mapper, see step 3 of `run_from_fastq.sh` |
| `--index` | yes | Bowtie2 prefix of the gene-tagged features (`mva_extract_features.py` + `bowtie2-build`) |
| `--large-n` | yes | Number of large cells (top barcodes by UMI rank) |
| `--small-n` | yes | Number of downsampled cells |
| `--empty-n` | yes | Number of simulated empty droplets |

`--large-n`, `--small-n` and `--empty-n` are all mandatory: they set the size of each
class and therefore the per-class recall that the summary reports.

## Usage

```bash
cd benchmarks
./run_benchmark.sh --sample 152 \
  --dedup-fq dedup.fq \
  --gff ea.gff3 --kegg-map gene_module.tsv --index bowtie2ref/mix3 \
  --large-n 1000 --small-n 2000 --empty-n 1000 \
  -t 32 -o results
```

Check the command sequence without running anything:

```bash
./run_benchmark.sh ... --dry-run
```

Reference files are cached in `-o`; rerunning the same sample skips steps 1–3 unless
`--force` is given. Collect several samples (skipping repeated headers):

```bash
awk 'FNR==1 && NR!=1 {next} {print}' results/*.summary.tsv | column -t -s $'\t'
```

### Options

| Option | Default | Meaning |
|---|---|---|
| `--reps INT` | 3 | Replicates per sample |
| `--empty-mean INT` | 100 | Mean reads per simulated empty barcode |
| `--min-ind INT` | 100000 | First rank index taken as ambient; must exceed `--large-n` |
| `--coef-min`, `--coef-max` | 0.2, 0.5 | Downsampling coefficient range |
| `--kneedle-s FLOAT` | 0.7 | Kneedle sensitivity `S` |
| `--pvalue FLOAT` | 0.01 | MVAcaller combined adjusted-P threshold |
| `--ind-min`, `--ind-max` | 300000, 865739999 | Ambient-profile rank window |
| `-t, --threads INT` | 8 | Aligner / samtools / MVAcaller threads |
| `--db FILE` | per replicate | Reuse one existing SQLite DB (it must have been built from the same mixed library); by default each replicate builds its own `molecules_<sample>_rep<N>.sqlite` |
| `--keep` | off | Keep mixed/simulated FASTQ and BAM files |
| `--force` | off | Rebuild the cached large/ambient reference files |
| `--dry-run` | off | Print commands only |

## Output

`results/<sample>.summary.tsv`, one row per method per replicate:

| Column | Meaning |
|---|---|
| `sample`, `rep` | Sample ID and replicate |
| `method` | `MVAcaller` or `knee` |
| `large`, `small`, `empty` | Barcodes of that class called as cells |
| `fdr` | `empty / (large + small + empty)`; the fraction of calls that are empties |

Reference files kept in `results/` (built once):

* `<sample>.large.fq`, `<sample>.large.labeled.fq`, `<sample>.large.sorted.bam`;
* `<sample>.ambient.fq`.

Per replicate:

* `test_<sample>_rep<N>.barcodes.tsv` — full MVAcaller table (per-barcode p-values,
  `call`, gene counts);
* `molecules_<sample>_rep<N>.sqlite` — molecule DB, reused if a replicate is re-run;
* with `--keep`, `empty_*.fq`, `small_*.fq`, `mix_*.fq`, `mix_*.sorted.bam`,
  `counts_*.tsv`.

## Differences from the manuscript script (`run_pipeline2.sh`)

The original script hard-coded `$HOME/py/*` lab utilities and fixed input paths; this version
parameterises them and uses the public tools shipped in `tools/`:

* the large/ambient split is now an explicit first step (`split_by_rank_index.py`), so the
  `rank_targ.py` groups and their `--target-umi-file` hand-off are gone;
* `simulate_empty2.py` and `gene_proportional_sampley.py` have no `--target-umi-file`:
  empty droplets are specified by count and mean reads (`--empty-n`, `--empty-mean`), and
  small cells by `--small-n` sampled with replacement inside the coefficient range;
* large/ambient are built once per sample instead of per replicate, matching how the
  original script reused fixed `large*.fq` / `ambient.fq` inputs;
* MVAcaller flags follow the release CLI (`--mode-omics` instead of `--mode omics_test`,
  `--build-sqlite` instead of `--sqlite-db`);
* summary statistics are taken by **column name** (`summarize_calls.awk`) instead of fixed
  column indices, so they survive table-format changes;
* each replicate builds its own SQLite DB; sharing one DB across replicates would score
  replicate 2/3 against replicate 1's molecules;
* `--qname-simple` is not passed because it is already the default; the inputs must carry
  BC/UMI in the read name (`..._<BC>_<UMI>`).

## Requirements

`python3` (numpy, pandas, scipy), `bowtie2`, `samtools`. `--dry-run` works without them.
