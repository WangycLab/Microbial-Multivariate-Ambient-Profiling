#!/usr/bin/env python3
"""
Build an MVAcaller molecules SQLite database from a count matrix instead of a BAM.

Use this when your reads were aligned/counted by another tool and you only have a
matrix (10x-style MTX, .h5ad, or a long-format table). The output has exactly the
schema produced by `MVAcaller.py --build-sqlite`:

    molecules(cb TEXT, umi TEXT, gene TEXT)

so the main program can consume it directly:

    python MVAcaller.py --build-sqlite sample.sqlite --reuse-sqlite \
        -i sample.sqlite -g annotation.gff3 --kegg-map gene_module.tsv -o out/sample

Performance
-----------
Molecules are not expanded one by one in Python. The matrix is staged into a
lightweight `cnt_stage(cb, gene, cnt)` table, and the expansion into one row per
molecule (`c1`, `c2`, ...) is done inside SQLite with a recursive CTE, one pass per
distinct count value.

Measured on a laptop SSD for an MTX with 1M non-zero entries / 2.5M molecules:
~4.5s for parsing plus expansion, ~5s to build the two indexes, ~10s in total
(~0.25M molecules/s end-to-end). The index build runs in C and does not depend on
the parser, so optimizing the Python side further cannot change it.

Deduplication
-------------
* Long-format tables that carry a UMI column are deduplicated on (cb, gene, umi),
  keeping the first occurrence - the same rule as `mva_dedup.py`.
* A pure count matrix (MTX / .h5ad) has no UMI information, so reads cannot be
  de-duplicated cell-wise. Each count unit becomes one molecule with a synthetic
  UMI. This is only valid if the counts are already UMI-collapsed molecules.

Requirements for the statistics to stay valid:
* The matrix must contain **empty droplets** if you rely on the automatic ambient
  profile (`--ind-min`). A cell-filtered matrix has no ambient signal to learn from.
* The matrix should retain **rRNA genes**, otherwise the rRNA fraction test degenerates.

Usage:
    python mva_mtx2sqlite.py --mtx matrix.mtx.gz --barcodes barcodes.tsv.gz \
        --features features.tsv.gz -o sample.sqlite
    python mva_mtx2sqlite.py --h5ad raw.h5ad -o sample.sqlite
    python mva_mtx2sqlite.py --table molecules.tsv -o sample.sqlite
"""

import argparse
import gzip
import os
import sqlite3
import sys
import time

SCHEMA = "CREATE TABLE molecules (cb TEXT NOT NULL, umi TEXT NOT NULL, gene TEXT NOT NULL)"
STAGING = "CREATE TABLE cnt_stage (cb TEXT, gene TEXT, cnt INTEGER)"
STAGING_INDEX = "CREATE INDEX idx_stage_cnt ON cnt_stage(cnt)"
INDEXES = (
    "CREATE INDEX idx_molecules_cb ON molecules(cb)",
    "CREATE INDEX idx_molecules_cb_gene ON molecules(cb, gene)",
)


def open_maybe_gzip(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def create_db(path):
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    conn = sqlite3.connect(path)
    cur = conn.cursor()
    cur.execute("PRAGMA journal_mode=WAL")
    cur.execute("PRAGMA synchronous=OFF")
    cur.execute("PRAGMA temp_store=MEMORY")
    cur.execute("PRAGMA cache_size=-200000")
    cur.execute("DROP TABLE IF EXISTS molecules")
    cur.execute(SCHEMA)
    cur.execute("DROP TABLE IF EXISTS cnt_stage")
    cur.execute(STAGING)
    cur.execute(STAGING_INDEX)
    conn.commit()
    return conn


def expand_staging(conn):
    """Turn one (cb, gene, cnt) row into `cnt` molecule rows, inside SQLite."""
    cur = conn.cursor()
    values = [count for (count,) in cur.execute("SELECT DISTINCT cnt FROM cnt_stage ORDER BY cnt")]
    for count in values:
        cur.execute(
            f"""
            WITH RECURSIVE nums(r) AS (SELECT 1 UNION ALL SELECT r+1 FROM nums WHERE r < {count})
            INSERT INTO molecules
            SELECT cb, 'c'||r, gene FROM cnt_stage, nums WHERE cnt = {count}
            """
        )
    conn.commit()


def finalize_db(conn):
    cur = conn.cursor()
    expand_staging(conn)
    cur.execute("DROP TABLE IF EXISTS cnt_stage")
    for statement in INDEXES:
        cur.execute(statement)
    conn.commit()
    conn.close()


def flush_stage(cur, batch):
    if batch:
        cur.executemany("INSERT INTO cnt_stage VALUES (?, ?, ?)", batch)
        batch.clear()


def flush_molecules(cur, batch):
    if batch:
        cur.executemany("INSERT INTO molecules VALUES (?, ?, ?)", batch)
        batch.clear()


def load_lines(path):
    with open_maybe_gzip(path) as handle:
        return [line.rstrip("\n").split("\t")[0] for line in handle if line.strip()]


# --------------------------------------------------------------------------- MTX
def read_mtx(conn, mtx_path, barcodes, features, strip_suffix, batch_size):
    """Stream a Matrix Market file (genes x cells, 1-based) into the staging table."""
    cur = conn.cursor()
    batch = []
    molecules = 0
    with open_maybe_gzip(mtx_path) as handle:
        dimensions_seen = False
        for line in handle:
            if line.startswith("%"):
                continue
            fields = line.split()
            if not dimensions_seen:
                dimensions_seen = True  # rows cols nnz
                continue
            if len(fields) < 3:
                continue
            gene_idx, cell_idx, value = int(fields[0]) - 1, int(fields[1]) - 1, float(fields[2])
            count = int(round(value))
            if count <= 0:
                continue
            cb = barcodes[cell_idx]
            if strip_suffix and cb.endswith(strip_suffix):
                cb = cb[:-len(strip_suffix)]
            batch.append((cb, features[gene_idx], count))
            molecules += count
            if len(batch) >= batch_size:
                flush_stage(cur, batch)
    flush_stage(cur, batch)
    conn.commit()
    return molecules


# -------------------------------------------------------------------------- H5AD
def read_h5ad(conn, h5ad_path, batch_size, chunk_cells, strip_suffix):
    """Stream an AnnData file cell-block by cell-block into the staging table."""
    try:
        import anndata as ad
        import numpy as np
        from scipy import sparse
    except ImportError:
        sys.exit("Reading .h5ad requires anndata, numpy and scipy: pip install anndata")

    adata = ad.read_h5ad(h5ad_path, backed="r")
    genes = list(adata.var_names)
    cells = list(adata.obs_names)
    cur = conn.cursor()
    batch = []
    molecules = 0
    n_cells = len(cells)

    for lo in range(0, n_cells, chunk_cells):
        hi = min(lo + chunk_cells, n_cells)
        block = adata.X[lo:hi]
        if sparse.issparse(block):
            block = block.tocsr()
            for row in range(block.shape[0]):
                cb = cells[lo + row]
                if strip_suffix and cb.endswith(strip_suffix):
                    cb = cb[:-len(strip_suffix)]
                start, end = block.indptr[row], block.indptr[row + 1]
                for col, value in zip(block.indices[start:end], block.data[start:end]):
                    count = int(round(float(value)))
                    if count <= 0:
                        continue
                    batch.append((cb, genes[col], count))
                    molecules += count
                if len(batch) >= batch_size:
                    flush_stage(cur, batch)
        else:
            arr = np.asarray(block)
            for row in range(arr.shape[0]):
                cb = cells[lo + row]
                if strip_suffix and cb.endswith(strip_suffix):
                    cb = cb[:-len(strip_suffix)]
                for col in np.nonzero(arr[row])[0]:
                    count = int(round(float(arr[row, col])))
                    if count <= 0:
                        continue
                    batch.append((cb, genes[col], count))
                    molecules += count
                if len(batch) >= batch_size:
                    flush_stage(cur, batch)

    flush_stage(cur, batch)
    conn.commit()
    return molecules, n_cells


# ------------------------------------------------------------------------- TABLE
def detect_columns(header):
    lower = [name.strip().lower() for name in header]

    def find(*candidates):
        for candidate in candidates:
            if candidate in lower:
                return lower.index(candidate)
        return None

    return {
        "cb": find("cb", "barcode", "cell", "cell_id"),
        "gene": find("gene", "gene_id", "feature", "geneid"),
        "umi": find("umi", "umi_id", "molecule"),
        "count": find("count", "counts", "n", "reads"),
    }


def read_table(conn, table_path, dedup, batch_size, sep, columns_spec):
    cur = conn.cursor()
    with open_maybe_gzip(table_path) as handle:
        header_line = handle.readline().rstrip("\n").lstrip("\ufeff")
        if sep is None:
            sep = "\t" if "\t" in header_line else ","
        header = [name.strip() for name in header_line.split(sep)]
        lower = [name.lower() for name in header]

        if columns_spec:
            names = [name.strip() for name in columns_spec.split(",")]
            try:
                cols = {"cb": lower.index(names[0].lower()),
                        "gene": lower.index(names[1].lower()),
                        "umi": None, "count": None}
                for name in names[2:]:
                    key = "umi" if name.lower().startswith("umi") else "count"
                    cols[key] = lower.index(name.lower())
            except (IndexError, ValueError):
                sys.exit(f"Column not found in header {header} for --columns {columns_spec}")
        else:
            cols = detect_columns(header)
            if cols["cb"] is None or cols["gene"] is None:
                sys.exit(
                    "Could not find barcode and gene columns. Header was: " + header_line
                    + "\nUse --columns cb,gene[,umi][,count] to name them explicitly."
                )

        seen = set()
        stage_batch = []
        molecule_batch = []
        rows = 0
        molecules = 0
        duplicates = 0
        max_index = max(index for index in cols.values() if index is not None)

        for line in handle:
            fields = line.rstrip("\n").split(sep)
            if len(fields) <= max_index:
                continue
            cb = fields[cols["cb"]]
            gene = fields[cols["gene"]]
            rows += 1

            if cols["umi"] is not None:
                umi = fields[cols["umi"]]
                if dedup:
                    key = (cb, gene, umi)
                    if key in seen:
                        duplicates += 1
                        continue
                    seen.add(key)
                molecule_batch.append((cb, umi, gene))
                molecules += 1
                if len(molecule_batch) >= batch_size:
                    flush_molecules(cur, molecule_batch)
            else:
                count = int(float(fields[cols["count"]])) if cols["count"] is not None else 1
                if count > 0:
                    stage_batch.append((cb, gene, count))
                    molecules += count
                if len(stage_batch) >= batch_size:
                    flush_stage(cur, stage_batch)

        flush_molecules(cur, molecule_batch)
        flush_stage(cur, stage_batch)
        conn.commit()

    has_umi = cols["umi"] is not None
    if has_umi and dedup:
        print(f"[dedup] {duplicates:,} duplicate (cb, gene, umi) rows removed", file=sys.stderr)
    return rows, molecules, cols


def main():
    parser = argparse.ArgumentParser(
        description="Build an MVAcaller molecules SQLite database from a count matrix "
                    "(MTX, .h5ad) or a long-format molecule table, without touching the BAM."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--mtx", help="Matrix Market file (10x style: rows = genes, cols = cells)")
    source.add_argument("--h5ad", help="AnnData .h5ad file")
    source.add_argument("--table", help="Long-format table with barcode/gene columns")

    parser.add_argument("--barcodes", help="barcodes.tsv[.gz] (required with --mtx)")
    parser.add_argument("--features", help="features.tsv[.gz] or genes.tsv[.gz] (required with --mtx)")
    parser.add_argument("--features-column", type=int, default=2, choices=(1, 2, 3),
                        help="Column of features.tsv to use as gene ID (default: 2; 3 = gene symbol)")
    parser.add_argument("--strip-suffix", default=None,
                        help="Remove this suffix from barcodes, e.g. '-1' for 10x barcodes")
    parser.add_argument("--columns", default=None,
                        help="Explicit column names for --table, e.g. 'cb,gene,umi' or "
                             "'cell,gene_id,count'")
    parser.add_argument("--sep", default=None, help="Field separator for --table (auto-detected)")
    parser.add_argument("--dedup", action=argparse.BooleanOptionalAction, default=True,
                        help="Deduplicate (cb, gene, umi) triplets when a UMI column exists "
                             "(default: on)")
    parser.add_argument("-o", "--output", default=None,
                        help="Output SQLite path (default: <input>.sqlite)")
    parser.add_argument("--batch-size", type=int, default=200000,
                        help="Rows per INSERT batch (default: 200000)")
    parser.add_argument("--chunk-cells", type=int, default=20000,
                        help="Cells per block when streaming .h5ad (default: 20000)")
    args = parser.parse_args()

    input_path = args.mtx or args.h5ad or args.table
    output = args.output or os.path.splitext(input_path)[0] + ".sqlite"

    started = time.time()
    conn = create_db(output)
    has_umi = False

    if args.mtx:
        if not args.barcodes or not args.features:
            sys.exit("--mtx requires --barcodes and --features")
        barcodes = load_lines(args.barcodes)
        feature_lines = [line.rstrip("\n").split("\t") for line in open_maybe_gzip(args.features)]
        column = args.features_column - 1
        features = [fields[column] if len(fields) > column else fields[1] for fields in feature_lines]
        n_molecules = read_mtx(conn, args.mtx, barcodes, features, args.strip_suffix,
                               args.batch_size)
        n_cells = len(barcodes)
        n_features = len(features)
    elif args.h5ad:
        n_molecules, n_cells = read_h5ad(conn, args.h5ad, args.batch_size,
                                         args.chunk_cells, args.strip_suffix)
        n_features = -1
    else:
        rows, n_molecules, cols = read_table(conn, args.table, args.dedup,
                                             args.batch_size, args.sep, args.columns)
        n_cells = rows
        n_features = -1
        has_umi = cols["umi"] is not None

    finalize_db(conn)

    print(f"Molecules written : {n_molecules:,}")
    print(f"Barcodes          : {n_cells:,}")
    if n_features > 0:
        print(f"Genes             : {n_features:,}")
    print(f"SQLite            : {output}")
    print(f"Elapsed           : {time.time() - started:.1f}s")

    if not has_umi:
        print("\n[dedup] This input carries no UMI column, so reads could NOT be "
              "de-duplicated: every count unit became one molecule. This is correct only "
              "if the matrix values are already UMI-collapsed molecules.", file=sys.stderr)
    print("\n[ambient] Remember that the automatic ambient profile needs empty droplets in "
          "the matrix; a cell-filtered matrix will silently give a biased profile. "
          "Check what --ind-min resolves to before trusting the calls.", file=sys.stderr)


if __name__ == "__main__":
    main()
