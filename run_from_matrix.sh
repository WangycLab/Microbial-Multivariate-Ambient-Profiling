#!/usr/bin/env bash
# run_from_matrix.sh — end-to-end MVAcaller run from a count matrix (MTX / .h5ad / table)
#
# Order of execution:
#   1. annotation -> GFF3 (GTF input is converted); rRNA genes are read from it
#   2. gene -> KEGG module map (use --kegg-map, or build it from --eggnog-db + --genome)
#   3. matrix -> SQLite molecule DB
#   4. MVAcaller cell calling (reuses the SQLite DB)
#
# No alignment is performed, so no Bowtie2 index is built. The matrix must still contain
# the ambient (empty-droplet) barcodes: MVAcaller learns the ambient profile from the
# low-UMI tail of the barcode rank curve, and a cell-filtered matrix silently produces a
# biased profile.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MTX=""; BARCODES=""; FEATURES=""; H5AD=""; TABLE=""
ANN=""; FNA=""; KM=""; OUT=""; EGGNOG_DB=""
STRIP=""; FEAT_COL=2
THREADS=8; MODE_OMICS=0; NO_DEDUP=0; DRY_RUN=0; FORCE=0

log()  { printf '\n[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die()  { printf '\nERROR: %s\n' "$*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "missing dependency: $1"; }

usage() {
cat <<'EOF'
Usage: run_from_matrix.sh --mtx matrix.mtx.gz --barcodes barcodes.tsv.gz \
                          --features features.tsv.gz -g annotation.gff3 \
                          -n genome.fna -o out/sample [options]

Input (choose one)
  --mtx FILE             Matrix Market file (10x raw_feature_bc_matrix)
  --barcodes FILE        barcodes.tsv(.gz) for --mtx
  --features FILE        features/genes.tsv(.gz) for --mtx
  --h5ad FILE            AnnData .h5ad with raw counts (needs anndata)
  --table FILE           Long table: barcode <tab> gene <tab> [umi] [count]

Required
  -g, --annotation FILE  GFF3 (or GTF) annotation
  -o, --out DIR          Output directory

Reference build
  -n, --genome FILE      Genome FASTA used to build the gene -> KEGG map (default path);
                         not needed when --kegg-map or a cached map is available
  --kegg-map FILE        Reuse an existing gene -> KEGG map instead of building one
  --eggnog-db DIR        eggNOG database directory for emapper.py
                         (default: the one configured for your eggNOG-mapper install)

Matrix details
  --strip-suffix TEXT    Strip this suffix from barcodes (e.g. "-1")
  --features-column INT  Column 1|2|3 of features.tsv (2=gene ID, 3=symbol)   [2]
  --no-dedup             Keep duplicate (barcode, gene, UMI) rows in --table input

Run control
  -t, --threads INT      Threads for SQLite build and statistics                 [8]
  --mode-omics           Multi-species (total-RNA third test); --umi-lower becomes 100
  --force                Recompute steps whose outputs already exist
  --dry-run              Print every command without executing it
  -h, --help             Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  # accept --opt=value as well as --opt value
  if [[ "$1" == --*=* ]]; then set -- "${1%%=*}" "${1#*=}" "${@:2}"; fi
  case "$1" in
    --mtx)              MTX="$2"; shift 2 ;;
    --barcodes)         BARCODES="$2"; shift 2 ;;
    --features)         FEATURES="$2"; shift 2 ;;
    --h5ad)             H5AD="$2"; shift 2 ;;
    --table)            TABLE="$2"; shift 2 ;;
    -g|--annotation)    ANN="$2"; shift 2 ;;
    -n|--genome)        FNA="$2"; shift 2 ;;
    --kegg-map)         KM="$2"; shift 2 ;;
    -o|--out)           OUT="$2"; shift 2 ;;
    --eggnog-db)        EGGNOG_DB="$2"; shift 2 ;;
    --strip-suffix)     STRIP="$2"; shift 2 ;;
    --features-column)  FEAT_COL="$2"; shift 2 ;;
    -t|--threads)       THREADS="$2"; shift 2 ;;
    --mode-omics)       MODE_OMICS=1; shift ;;
    --no-dedup)         NO_DEDUP=1; shift ;;
    --force)            FORCE=1; shift ;;
    --dry-run)          DRY_RUN=1; shift ;;
    -h|--help)          usage; exit 0 ;;
    *)                  die "unknown option: $1 (see --help)" ;;
  esac
done

[[ -n "$ANN" && -n "$OUT" ]] || { usage; die "missing required arguments"; }
[[ -n "$MTX" || -n "$H5AD" || -n "$TABLE" ]] || die "provide --mtx, --h5ad or --table"
[[ -z "$MTX" || ( -n "$BARCODES" && -n "$FEATURES" ) ]] \
  || die "--mtx also needs --barcodes and --features"
if [[ $DRY_RUN == 0 ]]; then need python3; fi

WORK="$OUT/work"
if [[ $FORCE == 1 && $DRY_RUN == 0 ]]; then rm -rf "$WORK"; fi
mkdir -p "$WORK"

# ---------------------------------------------------------------- 1. annotation -> GFF3
case "$ANN" in
  *.gtf|*.gtf.gz)
    log "Step 1 | GTF -> GFF3"
    GFF="$WORK/annotation.gff3"
    if [[ -s "$GFF" && $FORCE == 0 ]]; then
      log "reusing $GFF"
    else
      [[ $DRY_RUN == 1 ]] && log "python3 mva_gtf2gff3.py $ANN $GFF" \
        || python3 "$SCRIPT_DIR/mva_gtf2gff3.py" "$ANN" "$GFF"
    fi
    ;;
  *) GFF="$ANN" ;;
esac

# ------------------------------------------------------------- 2. gene -> KEGG map
log "Step 2 | gene -> KEGG module map"
if [[ -n "$KM" ]]; then
  [[ -r "$KM" ]] || die "cannot read: $KM"
  log "using $KM"
elif [[ -s "$WORK/gene_module.tsv" && $FORCE == 0 ]]; then
  KM="$WORK/gene_module.tsv"; log "reusing $KM"
else
  [[ -n "$FNA" ]] || die "building the KEGG map needs --genome (-n), or pass --kegg-map"
  [[ $DRY_RUN == 1 ]] || { need seqkit; need emapper.py; }
  KM="$WORK/gene_module.tsv"
  REF="$WORK/ref"
  EGGNOG_ARGS=(--itype proteins -o "$REF/eggnog" --cpu "$THREADS" --dbmem)
  if [[ -n "$EGGNOG_DB" ]]; then EGGNOG_ARGS+=(--data_dir "$EGGNOG_DB"); fi
  if [[ $DRY_RUN == 1 ]]; then
    log "python3 mva_extract_features.py -g $GFF -f $FNA -o $REF"
    log "seqkit translate $REF/cds.fa > $REF/protein.faa"
    log "emapper.py -i $REF/protein.faa ${EGGNOG_ARGS[*]}"
    log "grep -v '^#' $REF/eggnog.emapper.annotations | cut -f1,14 | sed 's/^[^#]*#//' > $KM"
  else
    python3 "$SCRIPT_DIR/mva_extract_features.py" -g "$GFF" -f "$FNA" -o "$REF"
    seqkit translate "$REF/cds.fa" > "$REF/protein.faa"
    emapper.py -i "$REF/protein.faa" "${EGGNOG_ARGS[@]}"
    grep -v '^#' "$REF/eggnog.emapper.annotations" | cut -f1,14 \
      | sed 's/^[^#]*#//' > "$KM"
  fi
fi

# --------------------------------------------------- 3. matrix -> SQLite
log "Step 3 | count matrix -> SQLite molecule DB"
SQLITE="$WORK/molecules.sqlite"
if [[ -s "$SQLITE" && $FORCE == 0 ]]; then
  log "reusing $SQLITE"
else
  CONV_ARGS=(-o "$SQLITE")
  if [[ -n "$MTX" ]]; then
    CONV_ARGS+=(--mtx "$MTX" --barcodes "$BARCODES" --features "$FEATURES"
                --features-column "$FEAT_COL")
  elif [[ -n "$H5AD" ]]; then
    CONV_ARGS+=(--h5ad "$H5AD")
  else
    CONV_ARGS+=(--table "$TABLE")
    [[ $NO_DEDUP == 1 ]] && CONV_ARGS+=(--no-dedup)
  fi
  [[ -n "$STRIP" ]] && CONV_ARGS+=(--strip-suffix="$STRIP")
  log "python3 mva_mtx2sqlite.py ${CONV_ARGS[*]}"
  [[ $DRY_RUN == 1 ]] || python3 "$SCRIPT_DIR/mva_mtx2sqlite.py" "${CONV_ARGS[@]}"
fi

# ------------------------------------------------------- 4. cell calling
log "Step 4 | cell calling"
MVA_ARGS=(-g "$GFF" --kegg-map "$KM" --build-sqlite "$SQLITE" --reuse-sqlite
          -o "$OUT/mva_run" --threads "$THREADS")
if [[ $MODE_OMICS == 1 ]]; then MVA_ARGS+=(--mode-omics); fi
log "python3 MVAcaller.py ${MVA_ARGS[*]}"
[[ $DRY_RUN == 1 ]] || python3 "$SCRIPT_DIR/MVAcaller.py" "${MVA_ARGS[@]}"

log "Done."
[[ $DRY_RUN == 1 ]] && exit 0
printf 'cells    : %s\n' "$OUT/mva_run.cells.txt"
printf 'table    : %s\n' "$OUT/mva_run.barcodes.tsv"
printf 'summary  : %s\n' "$OUT/mva_run.summary.txt"
printf 'n cells  : %s\n' "$(wc -l < "$OUT/mva_run.cells.txt")"
