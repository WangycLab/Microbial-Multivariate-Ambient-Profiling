#!/usr/bin/env python3
"""
split_by_rank_index.py
======================
按 barcode UMI count 降序排序，根据 max-ind / min-ind 将 reads 分成 large / ambient。

用法:
  # 方式一：从头扫描 FASTQ（按 read count 降序排序）
  python split_by_rank_index.py -i dedup.fq --max-ind 500 -l large.fq

  # 方式二：同时指定 max-ind 和 min-ind
  python split_by_rank_index.py -i dedup.fq --max-ind 500 --min-ind 2000 \
                                -l large.fq -a ambient.fq

  # 方式三：输入已有 counts.tsv（已按 umi_count 降序排列），跳过扫描排序
  python split_by_rank_index.py -i dedup.fq --counts counts.tsv --max-ind 500 \
                                -l large.fq

输出:
  large.fq   — UMI rank 1 ~ max_ind 的 barcode reads
  ambient.fq — UMI rank > min_ind 的 barcode reads（仅当指定 --min-ind 时输出）
"""

import sys
import gzip
import argparse
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(
        description='Split FASTQ by barcode UMI rank',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # Only --max-ind: output large.fq only\n"
            "  python split_by_rank_index.py -i dedup.fq --max-ind 500 -l large.fq\n\n"
            "  # Both --max-ind and --min-ind: output large + ambient\n"
            "  python split_by_rank_index.py -i dedup.fq --max-ind 500 --min-ind 2000 \\\n"
            "      -l large.fq -a ambient.fq\n\n"
            "  # Use pre-computed counts.tsv (skip scan & sort)\n"
            "  python split_by_rank_index.py -i dedup.fq --counts counts.tsv --max-ind 500 \\\n"
            "      -l large.fq\n"
        ),
    )
    p.add_argument('-i', '--input', required=True, help='Input FASTQ (.gz supported)')
    p.add_argument('--max-ind', type=int, default=None,
                   help='Max index: barcodes with rank <= max_ind go to large')
    p.add_argument('--min-ind', type=int, default=None,
                   help='Min index: barcodes with rank > min_ind go to ambient')
    p.add_argument('--counts', default=None,
                   help='Pre-computed counts.tsv (umi_count descending), skip first-pass scan')
    p.add_argument('-l', '--large', default='large.fq', help='Output large.fq')
    p.add_argument('-a', '--ambient', default='ambient.fq', help='Output ambient.fq (only with --min-ind)')
    return p.parse_args()


def load_counts(counts_path):
    """从 counts.tsv 加载已排好序的 barcode 列表。格式: barcode reads_count umi_count"""
    bc_sorted = []
    with open(counts_path) as f:
        header = f.readline()
        if not header.startswith('barcode'):
            print(f'[ERROR] counts.tsv header expected to start with "barcode", got: {header.rstrip()}',
                  file=sys.stderr)
            sys.exit(1)
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            bc, reads, umis = parts[0], int(parts[1]), int(parts[2])
            bc_sorted.append((bc, reads, umis))
    return bc_sorted


def scan_fastaq(input_path):
    """第一遍扫描 FASTQ，统计 barcode reads 数，按 read count 降序排序。"""
    opener = gzip.open if input_path.endswith('.gz') else open
    bc_counts = defaultdict(int)
    total_reads = 0

    with opener(input_path, 'rt') as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            total_reads += 1

            if not header.startswith('@'):
                continue
            parts = header.rstrip('\n').split('_')
            if len(parts) < 3:
                continue
            bc = parts[-2]
            bc_counts[bc] += 1

    bc_sorted = sorted(bc_counts.items(), key=lambda x: x[1], reverse=True)
    return bc_sorted, total_reads


def main():
    args = parse_args()

    if args.max_ind is None:
        print('[ERROR] --max-ind is required', file=sys.stderr)
        sys.exit(1)

    # === 获取 barcode 排序列表 ===
    if args.counts:
        print('[INFO] Loading barcode ranks from counts.tsv...', file=sys.stderr)
        bc_sorted = load_counts(args.counts)
        print(f'  Loaded {len(bc_sorted):,} barcodes (sorted by umi_count desc)',
              file=sys.stderr)
    else:
        print('[1/2] Counting reads per barcode...', file=sys.stderr)
        bc_sorted, total_reads = scan_fastaq(args.input)
        print(f'  Total reads: {total_reads:,}, Unique barcodes: {len(bc_sorted):,}',
              file=sys.stderr)

    if bc_sorted:
        # bc_sorted 可能是 (bc, count) 或 (bc, reads, umis)
        max_val = bc_sorted[0][1]
        min_val = bc_sorted[-1][1]
        print(f'  Count range: max={max_val:,}, min={min_val:,}', file=sys.stderr)

    # === 分类 barcode ===
    large_bcs = set()
    ambient_bcs = set()

    for rank, item in enumerate(bc_sorted, start=1):
        bc = item[0]
        if rank <= args.max_ind:
            large_bcs.add(bc)
        elif args.min_ind is not None and rank > args.min_ind:
            ambient_bcs.add(bc)
        # rank 在 max_ind+1 ~ min_ind 之间的 reads 不输出

    print(f'  Large:   {len(large_bcs):,} barcodes (rank 1~{args.max_ind})',
          file=sys.stderr)
    if args.min_ind is not None:
        print(f'  Ambient: {len(ambient_bcs):,} barcodes (rank > {args.min_ind})',
              file=sys.stderr)

    # === 第二遍：输出 ===
    print('\n[2/2] Writing output files...', file=sys.stderr)
    opener = gzip.open if args.input.endswith('.gz') else open

    fl = open(args.large, 'w')
    fa = None
    if args.min_ind is not None:
        if args.ambient:
            fa = open(args.ambient, 'w')
        else:
            print('[WARN] --min-ind set but no --ambient output file, ambient reads will be discarded',
                  file=sys.stderr)

    counts_large = 0
    counts_ambient = 0
    counts_skipped = 0

    with opener(args.input, 'rt') as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()

            if not header.startswith('@'):
                continue
            parts = header.rstrip('\n').split('_')
            if len(parts) < 3:
                continue
            bc = parts[-2]

            if bc in large_bcs:
                fl.write(header)
                fl.write(seq)
                fl.write(plus)
                fl.write(qual)
                counts_large += 1
            elif fa is not None and bc in ambient_bcs:
                fa.write(header)
                fa.write(seq)
                fa.write(plus)
                fa.write(qual)
                counts_ambient += 1
            else:
                counts_skipped += 1

    fl.close()
    if fa:
        fa.close()

    print(f'  {args.large}:   {counts_large:,} reads', file=sys.stderr)
    if fa:
        print(f'  {args.ambient}: {counts_ambient:,} reads', file=sys.stderr)
    if counts_skipped:
        print(f'  Skipped (mix region): {counts_skipped:,} reads', file=sys.stderr)
    print('\nDone!', file=sys.stderr)


if __name__ == '__main__':
    main()
