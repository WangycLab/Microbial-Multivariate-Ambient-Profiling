#!/usr/bin/env python3
"""
Load sparse TSV (gene\tcell\tcount) into scanpy-compatible AnnData h5ad.
Two-pass chunked reading, never materializes a dense matrix.

Output: obs=cells, var=genes, X=csr_matrix (cells × genes) — matches scanpy convention.

Usage:
  python mva_load_h5ad.py -i wt.knee.tsv -o wt.h5ad
"""

import argparse
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad


def main():
    parser = argparse.ArgumentParser(description="Load sparse TSV to scanpy-compatible h5ad")
    parser.add_argument("-i", "--input", required=True, help="Sparse TSV: gene\\tcell\\tcount")
    parser.add_argument("-o", "--output", required=True, help="Output .h5ad")
    args = parser.parse_args()

    chunk_size = 500000

    # ── Pass 1: collect unique gene & cell names ──
    print("[1/3] Scanning genes and cells...", file=sys.stderr, flush=True)
    gene_set = set()
    cell_set = set()
    with open(args.input) as f:
        next(f)  # header
        for line in f:
            parts = line.split("\t", 2)
            if len(parts) < 3:
                continue
            gene_set.add(parts[0])
            cell_set.add(parts[1].rstrip("\n"))

    genes = sorted(gene_set)
    cells = sorted(cell_set)
    gene_idx = {g: i for i, g in enumerate(genes)}
    cell_idx = {c: i for i, c in enumerate(cells)}
    print(f"  Genes: {len(genes):,}  Cells: {len(cells):,}", file=sys.stderr, flush=True)

    # ── Pass 2: read into COO arrays via chunks ──
    print("[2/3] Reading counts...", file=sys.stderr, flush=True)
    row_g = []   # gene index
    col_c = []   # cell index
    data_v = []  # count
    total = 0

    for chunk in pd.read_csv(
        args.input, sep="\t", header=0,
        names=["gene", "cell", "count"],
        dtype={"gene": "category", "cell": "category", "count": np.int32},
        chunksize=chunk_size,
    ):
        if chunk.empty:
            continue
        g_arr = chunk["gene"].astype(str).values
        c_arr = chunk["cell"].astype(str).values
        v_arr = chunk["count"].values
        for g, c, v in zip(g_arr, c_arr, v_arr):
            gi = gene_idx.get(g)
            ci = cell_idx.get(c)
            if gi is not None and ci is not None:
                row_g.append(gi)
                col_c.append(ci)
                data_v.append(v)
        total += len(chunk)
        if total % (chunk_size * 10) == 0:
            print(f"  Read {total:,} rows...", file=sys.stderr, flush=True)

    print(f"  Non-zero entries: {len(data_v):,}", file=sys.stderr, flush=True)

    # ── Build csr: gene × cell, then transpose to cell × gene ──
    print("[3/3] Building AnnData...", file=sys.stderr, flush=True)
    X_gc = sp.csr_matrix(
        (np.array(data_v, dtype=np.int32),
         (np.array(row_g, dtype=np.int32),
          np.array(col_c, dtype=np.int32))),
        shape=(len(genes), len(cells))
    )
    # Transpose: cell × gene, convert to csr for scanpy
    X = sp.csr_matrix(X_gc.T)

    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=cells),
        var=pd.DataFrame(index=genes),
    )
    adata.obs_names = cells
    adata.var_names = genes
    adata.layers["counts"] = X.copy()

    adata.write_h5ad(args.output)
    print(f"  Saved: {args.output}", file=sys.stderr, flush=True)
    print(f"  Shape: {X.shape} (cells × genes), nnz={X.nnz:,}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
