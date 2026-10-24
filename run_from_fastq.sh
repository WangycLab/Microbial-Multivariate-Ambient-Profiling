#!/usr/bin/env bash
# run_from_fastq.sh — end-to-end MVAcaller run from paired FASTQ (R1/R2) + genome/annotation
#
# Order of execution:
#   1. annotation -> GFF3 (GTF input is converted)
#   2. GFF3 + genome -> gene-tagged FASTA (CDS/tRNA/rRNA/ncRNA) and Bowtie2 index
#   3. gene -> KEGG module map (use --kegg-map, or build it from --eggnog-db)
#   4. R1 BC/UMI -> R2 read header (skipped with --skip-extract)
#   5. align to the gene-tagged index -> coordinate-sorted BAM
#   6. UMI deduplication -> dedup.bam (+ counts.tsv)
#   7. SQLite molecule DB -> MVAcaller cell calling
#
# Only steps 5-7 are platform-independent. Step 4 assumes the barcode/UMI sit at fixed
# positions in R1; adjust --bc-start/--bc-len/--umi-start/--umi-len, or use --skip-extract
# and preprocess the reads yourself (MVAcaller then reads BC/UMI from the read name).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R1=""; R2=""; FNA=""; ANN=""; KM=""; OUT=""; EGGNOG_DB=""
BC_START=1; BC_LEN=20; UMI_START=21; UMI_LEN=8
THREADS=8; SKIP_EXTRACT=0; MODE_OMICS=0; DRY_RUN=0; FORCE=0; DENSE_MATRIX=0
MAT_OUT=""; MAT_FLAG=""

log()  { printf '\n[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die()  { printf '\nERROR: %s\n' "$*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "missing dependency: $1"; }

usage() {
cat <<'EOF'
Usage: run_from_fastq.sh -1 R1.fq.gz -2 R2.fq.gz -n genome.fna -g annotation.gff3 \
                         --kegg-map gene_module.tsv -o out/sample [options]

Required
  -1, --r1 FILE          R1 FASTQ (.gz ok) containing barcode + UMI
  -2, --r2 FILE          R2 FASTQ (.gz ok) containing the transcript
  -n, --genome FILE      Genome FASTA (.fna) from which features are extracted
  -g, --annotation FILE  GFF3 (or GTF) annotation for the same assembly
  -o, --out DIR          Output directory

Barcode / UMI layout (platform specific)
  --bc-start INT        1-based barcode start in R1            [1]
  --bc-len INT           Barcode length                         [20]
  --umi-start INT       1-based UMI start in R1                 [21]
  --umi-len INT          UMI length                             [8]
  --skip-extract         Skip step 4; the FASTQ read names already carry BC/UMI

Reference build
  --kegg-map FILE        Reuse an existing gene -> KEGG map instead of building one
  --eggnog-db DIR        eggNOG database directory for emapper.py
                         (default: the one configured for your eggNOG-mapper install)

Run control
  -t, --threads INT      Threads for aligner / dedup / SQLite / statistics   [8]
  --mode-omics           Multi-species (total-RNA third test); --umi-lower becomes 100
  --dense-matrix         Write a dense gene x cell table instead of sparse COO
  --force                Recompute steps whose outputs already exist
  --dry-run              Print every command without executing it
  -h, --help             Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  # accept --opt=value as well as --opt value
  if [[ "$1" == --*=* ]]; then set -- "${1%%=*}" "${1#*=}" "${@:2}"; fi
  case "$1" in
    -1|--r1)        R1="$2"; shift 2 ;;
    -2|--r2)        R2="$2"; shift 2 ;;
    -n|--genome)    FNA="$2"; shift 2 ;;
    -g|--annotation) ANN="$2"; shift 2 ;;
    --kegg-map)     KM="$2"; shift 2 ;;
    -o|--out)       OUT="$2"; shift 2 ;;
    --eggnog-db)    EGGNOG_DB="$2"; shift 2 ;;
    --bc-start)     BC_START="$2"; shift 2 ;;
    --bc-len)       BC_LEN="$2"; shift 2 ;;
    --umi-start)    UMI_START="$2"; shift 2 ;;
    --umi-len)      UMI_LEN="$2"; shift 2 ;;
    -t|--threads)   THREADS="$2"; shift 2 ;;
    --skip-extract) SKIP_EXTRACT=1; shift ;;
    --mode-omics)   MODE_OMICS=1; shift ;;
    --dense-matrix) DENSE_MATRIX=1; shift ;;
    --force)        FORCE=1; shift ;;
    --dry-run)      DRY_RUN=1; shift ;;
    -h|--help)      usage; exit 0 ;;
    *)              die "unknown option: $1 (see --help)" ;;
  esac
done

[[ -n "$R1" && -n "$R2" && -n "$FNA" && -n "$ANN" && -n "$OUT" ]] \
  || { usage; die "missing required arguments"; }
if [[ $SKIP_EXTRACT == 0 && $DRY_RUN == 0 ]]; then need go; fi
for f in "$R1" "$R2" "$FNA" "$ANN"; do
  [[ -r "$f" ]] || die "cannot read: $f"
done
if [[ $DRY_RUN == 0 ]]; then
  need bowtie2; need samtools; need python3
fi

WORK="$OUT/work"; REF="$OUT/ref"
if [[ $FORCE == 1 && $DRY_RUN == 0 ]]; then rm -rf "$WORK" "$REF"; fi
mkdir -p "$WORK" "$REF"

# ---------------------------------------------------------------- 1. annotation -> GFF3
case "$ANN" in
  *.gtf|*.gtf.gz)
    log "Step 1 | GTF -> GFF3"
    GFF="$REF/annotation.gff3"
    if [[ -s "$GFF" && $FORCE == 0 ]]; then
      log "reusing $GFF"
    else
      [[ $DRY_RUN == 1 ]] && log "python3 mva_gtf2gff3.py $ANN $GFF" \
        || python3 "$SCRIPT_DIR/mva_gtf2gff3.py" "$ANN" "$GFF"
    fi
    ;;
  *) GFF="$ANN" ;;
esac

# ------------------------------------------- 2. gene-tagged features + Bowtie2 index
log "Step 2 | gene-tagged features and Bowtie2 index"
if [[ -s "$REF/features.fa" && $FORCE == 0 ]]; then
  log "reusing $REF/features.fa"
else
  [[ $DRY_RUN == 1 ]] && log "python3 mva_extract_features.py -g $GFF -f $FNA -o $REF" \
    || python3 "$SCRIPT_DIR/mva_extract_features.py" -g "$GFF" -f "$FNA" -o "$REF"
fi
if [[ -s "$REF/features.1.bt2" && $FORCE == 0 ]]; then
  log "reusing $REF/features.*.bt2"
else
  [[ $DRY_RUN == 1 ]] && log "bowtie2-build --threads $THREADS $REF/features.fa $REF/features" \
    || bowtie2-build --threads "$THREADS" "$REF/features.fa" "$REF/features" >/dev/null
fi

# ------------------------------------------------------------- 3. gene -> KEGG map
log "Step 3 | gene -> KEGG module map"
if [[ -n "$KM" ]]; then
  [[ -r "$KM" ]] || die "cannot read: $KM"
  log "using $KM"
elif [[ -s "$REF/gene_module.tsv" && $FORCE == 0 ]]; then
  KM="$REF/gene_module.tsv"; log "reusing $KM"
else
  [[ $DRY_RUN == 1 ]] || { need seqkit; need emapper.py; }
  KM="$REF/gene_module.tsv"
  EGGNOG_ARGS=(--itype proteins -o "$REF/eggnog" --cpu "$THREADS" --dbmem)
  if [[ -n "$EGGNOG_DB" ]]; then EGGNOG_ARGS+=(--data_dir "$EGGNOG_DB"); fi
  if [[ $DRY_RUN == 1 ]]; then
    log "seqkit translate $REF/cds.fa > $REF/protein.faa"
    log "emapper.py -i $REF/protein.faa ${EGGNOG_ARGS[*]}"
    log "grep -v '^#' $REF/eggnog.emapper.annotations | cut -f1,14 | sed 's/^[^#]*#//' > $KM"
  else
    seqkit translate "$REF/cds.fa" > "$REF/protein.faa"
    emapper.py -i "$REF/protein.faa" "${EGGNOG_ARGS[@]}"
    grep -v '^#' "$REF/eggnog.emapper.annotations" | cut -f1,14 \
      | sed 's/^[^#]*#//' > "$KM"
  fi
fi

# ----------------------------------------------- 4. barcode / UMI -> R2 read header
log "Step 4 | barcode and UMI extraction"
if [[ $SKIP_EXTRACT == 1 ]]; then
  R2_USE="$R2"; log "skipped (--skip-extract); read names must be ..._<BC>_<UMI>"
else
  EXTRACT_BIN="$WORK/mva_umi_extract"
  [[ -x "$EXTRACT_BIN" ]] || {
    [[ $DRY_RUN == 1 ]] && log "go build -o $EXTRACT_BIN $SCRIPT_DIR/mva_umi_extract.go" \
      || go build -o "$EXTRACT_BIN" "$SCRIPT_DIR/mva_umi_extract.go"
  }
  R2_USE="$WORK/R2.extracted.fq.gz"
  [[ $DRY_RUN == 1 ]] && log "$EXTRACT_BIN -1 $R1 -2 $R2 -o $R2_USE --bc-start $BC_START --bc-len $BC_LEN --umi-start $UMI_START --umi-len $UMI_LEN -t $THREADS" \
    || "$EXTRACT_BIN" -1 "$R1" -2 "$R2" -o "$R2_USE" \
         --bc-start "$BC_START" --bc-len "$BC_LEN" \
         --umi-start "$UMI_START" --umi-len "$UMI_LEN" -t "$THREADS"
fi

# ------------------------------------------------------------- 5. align + sort
log "Step 5 | alignment to the gene-tagged index"
SORTED="$WORK/temp.sorted.bam"
if [[ -s "$SORTED" && $FORCE == 0 ]]; then
  log "reusing $SORTED"
elif [[ $DRY_RUN == 1 ]]; then
  log "bowtie2 -x $REF/features -U $R2_USE --very-sensitive-local -p $THREADS | samtools view -bS - | samtools sort -@ $THREADS -o $SORTED"
else
  bowtie2 -x "$REF/features" -U "$R2_USE" --very-sensitive-local -p "$THREADS" \
    | samtools view -bS - \
    | samtools sort -@ "$THREADS" -o "$SORTED" -
fi

# --------------------------------------------------------------- 6. UMI dedup
log "Step 6 | UMI deduplication"
DEDUP="$WORK/dedup.bam"
if [[ -s "$DEDUP" && $FORCE == 0 ]]; then
  log "reusing $DEDUP"
elif [[ $DRY_RUN == 1 ]]; then
  log "python3 mva_dedup.py -b $SORTED -c $WORK/counts.tsv -d $DEDUP -o $WORK/dedup.fq"
else
  python3 "$SCRIPT_DIR/mva_dedup.py" -b "$SORTED" -c "$WORK/counts.tsv" \
    -d "$DEDUP" -o "$WORK/dedup.fq"
fi

# -------------------------------------------- 7. SQLite + cell calling
log "Step 7 | SQLite molecule DB and cell calling"
SQLITE="$WORK/molecules.sqlite"
if [[ ! -s "$SQLITE" || $FORCE == 1 ]]; then
  if command -v go >/dev/null 2>&1; then
    BUILDER="$WORK/mva_build_sqlite"
    [[ -x "$BUILDER" ]] || {
      [[ $DRY_RUN == 1 ]] && log "go build -o $BUILDER $SCRIPT_DIR/mva_build_sqlite.go" \
        || go build -o "$BUILDER" "$SCRIPT_DIR/mva_build_sqlite.go"
    }
    [[ $DRY_RUN == 1 ]] && log "$BUILDER -i $DEDUP -o $SQLITE -t $THREADS" \
      || "$BUILDER" -i "$DEDUP" -o "$SQLITE" -t "$THREADS"
  else
    log "Go not found; MVAcaller will build the SQLite DB with its Python builder"
  fi
fi

MVA_ARGS=(-g "$GFF" --kegg-map "$KM" --build-sqlite "$SQLITE" --reuse-sqlite
          -o "$OUT/mva_run" --threads "$THREADS")
if [[ $MODE_OMICS == 1 ]]; then MVA_ARGS+=(--mode-omics); fi
if [[ -s "$SQLITE" ]]; then
  log "python3 MVAcaller.py ${MVA_ARGS[*]}"
  [[ $DRY_RUN == 1 ]] || python3 "$SCRIPT_DIR/MVAcaller.py" "${MVA_ARGS[@]}"
else
  log "python3 MVAcaller.py -i $DEDUP ${MVA_ARGS[*]}"
  [[ $DRY_RUN == 1 ]] || python3 "$SCRIPT_DIR/MVAcaller.py" -i "$DEDUP" "${MVA_ARGS[@]}"
fi

# ------------------------------- 8. filtered gene x cell matrix
log "Step 8 | filtered gene x cell matrix"
if [[ -s "$OUT/mva_run.cells.txt" ]]; then
  MAT="$OUT/mva_run.filtered_matrix"
  if [[ $DENSE_MATRIX == 1 ]]; then
    MAT_OUT="$MAT.tsv"; MAT_FLAG=""
  else
    MAT_OUT="$MAT.sparse.tsv"; MAT_FLAG="--sparse"
  fi
  log "python3 mva_bam_to_mat.py -i $DEDUP -b $OUT/mva_run.cells.txt -o $MAT_OUT -t $THREADS ${MAT_FLAG}"
  [[ $DRY_RUN == 1 ]] || python3 "$SCRIPT_DIR/mva_bam_to_mat.py" -i "$DEDUP" \
      -b "$OUT/mva_run.cells.txt" -o "$MAT_OUT" -t "$THREADS" ${MAT_FLAG}
else
  log "no called cells; skipping the matrix step"
fi

log "Done."
[[ $DRY_RUN == 1 ]] && exit 0
printf 'cells    : %s\n' "$OUT/mva_run.cells.txt"
printf 'table    : %s\n' "$OUT/mva_run.barcodes.tsv"
printf 'summary  : %s\n' "$OUT/mva_run.summary.txt"
printf 'matrix   : %s\n' "$MAT_OUT"
printf 'n cells  : %s\n' "$(wc -l < "$OUT/mva_run.cells.txt")"
