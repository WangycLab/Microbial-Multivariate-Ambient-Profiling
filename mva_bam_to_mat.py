#!/usr/bin/env python3

import argparse
import subprocess
import pandas as pd
from collections import defaultdict

parser = argparse.ArgumentParser(
    description="Build gene x cell matrix from deduplicated BAM (samtools version)."
)
parser.add_argument("-i", "--input", required=True,
                    help="Input deduplicated BAM")
parser.add_argument("-b", "--barcodes",
                    help="Cell barcode list (optional; if not provided, keep all barcodes)")
parser.add_argument("-o", "--output", required=True,
                    help="Output matrix file (.tsv dense, or .sparse.tsv COO)")
parser.add_argument("-t", "--threads", type=int, default=4,
                    help="samtools threads (default: 4)")
parser.add_argument("--sparse", action="store_true",
                    help="Output sparse COO format (gene, cell, count; non-zero only)")
args = parser.parse_args()

# 读取保留的barcode（如果提供了）
keep = None
if args.barcodes:
    with open(args.barcodes) as f:
        keep = set(line.strip() for line in f if line.strip())

counts = defaultdict(int)

# samtools view -F 4: 只输出mapped reads
# cut -f1,3:       只取 QNAME 和 RNAME 两列，减少 Python 解析量
cmd = (
    f"samtools view -F 4 --no-header -@ {args.threads} "
    f"'{args.input}' | cut -f1,3"
)

with subprocess.Popen(
    cmd, stdout=subprocess.PIPE, text=True, shell=True,
    bufsize=8 * 1024 * 1024,
) as proc:

    for line in proc.stdout:
        tab = line.index("\t")
        qname = line[:tab]
        ref = line[tab + 1:].rstrip("\n")

        if "#" not in ref:
            continue

        gene = ref.split("#", 1)[1]
        if not gene:
            continue

        parts = qname.split("_")
        if len(parts) < 3:
            continue

        cb = parts[-2]
        if keep is not None and cb not in keep:
            continue

        counts[(gene, cb)] += 1

    proc.wait()

if args.sparse:
    # Sparse COO format: gene, cell, count (non-zero entries only)
    with open(args.output, "w") as fout:
        fout.write("gene\tcell\tcount\n")
        for (gene, cb), cnt in sorted(counts.items()):
            fout.write(f"{gene}\t{cb}\t{cnt}\n")
    n_genes = len({g for g, _ in counts})
    n_cells = len({c for _, c in counts})
    print(f"[bam2mat] sparse: {n_genes:,} genes x {n_cells:,} cells, "
          f"{len(counts):,} non-zero entries -> {args.output}")
else:
    # Dense matrix
    df = pd.Series(counts).unstack(fill_value=0)

    if keep is not None:
        df = df.reindex(columns=sorted(keep), fill_value=0)

    df = df.sort_index()
    df.to_csv(args.output, sep="\t")
    print(f"[bam2mat] dense: {df.shape[0]:,} genes x {df.shape[1]:,} cells -> {args.output}")
