#!/usr/bin/env python3
"""
simulate_empty.py
=================
从 empty FASTQ 中随机抽取 reads，组合成带有 artificial barcodes 的 simulated FASTQ。

用法:
  python simulate_empty.py -i empty.fq.gz -o simulated.fq.gz \
                           -n 1000 -m 50 -v 100 --seed 42
"""

import argparse
import random
import sys
import gzip
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description='Simulate empty droplets from empty FASTQ')
    p.add_argument('-i', '--input', required=True, help='Input empty FASTQ (.fq or .fq.gz)')
    p.add_argument('-o', '--output', required=True, help='Output simulated FASTQ (.fq or .fq.gz)')
    p.add_argument('-n', '--n-barcodes', type=int, default=1000,
                   help='Number of artificial barcodes. Default: 1000')
    p.add_argument('-m', '--mean', type=float, default=100,
                   help='Mean reads per barcode. Default: 100')
    p.add_argument('-v', '--variance', type=float, default=None,
                   help='Variance of reads per barcode. Default: mean (Poisson-like)')
    p.add_argument('-d', '--distribution', choices=['poisson', 'nb', 'normal'],
                   default='nb', help='Distribution of reads per barcode. Default: nb')
    p.add_argument('--barcode-prefix', default='EMPTY',
                   help='Prefix for artificial barcodes. Default: EMPTY')
    p.add_argument('--seed', type=int, default=42)
    return p.parse_args()


def read_fastq(fq_path):
    """读取 FASTQ，返回 [(name, seq, plus, qual), ...]"""
    opener = gzip.open if fq_path.endswith('.gz') else open
    records = []
    with opener(fq_path, 'rt') as f:
        while True:
            name = f.readline().strip()
            if not name:
                break
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            records.append((name, seq, plus, qual))
    return records


def generate_barcode_sizes(n_barcodes, mean, variance, dist, seed):
    """生成每个 barcode 的 reads 数"""
    np.random.seed(seed)
    if variance is None:
        variance = mean

    if dist == 'poisson':
        sizes = np.random.poisson(mean, n_barcodes)
    elif dist == 'normal':
        std = np.sqrt(variance)
        sizes = np.random.normal(mean, std, n_barcodes)
        sizes = np.maximum(sizes, 1).astype(int)
    elif dist == 'nb':
        # 负二项分布: mean = r*(1-p)/p, var = r*(1-p)/p^2
        # -> r = mean^2 / (var - mean), p = r / (r + mean)
        if variance <= mean:
            # 方差 <= 均值时退化为泊松
            sizes = np.random.poisson(mean, n_barcodes)
        else:
            r = mean ** 2 / (variance - mean)
            p = r / (r + mean)
            sizes = np.random.negative_binomial(r, p, n_barcodes)
    else:
        raise ValueError(f"Unknown distribution: {dist}")

    return sizes.astype(int)


def main():
    args = parse_args()
    random.seed(args.seed)
    np.random.seed(args.seed)

    # 读取所有 reads
    print(f"Reading {args.input}...", file=sys.stderr)
    records = read_fastq(args.input)
    total_reads = len(records)
    print(f"  Total reads: {total_reads:,}", file=sys.stderr)

    if total_reads == 0:
        print("Error: no reads found.", file=sys.stderr)
        sys.exit(1)

    # 生成每个 barcode 的 reads 数
    sizes = generate_barcode_sizes(
        args.n_barcodes, args.mean, args.variance,
        args.distribution, args.seed
    )
    requested = int(sizes.sum())

    # 如果请求 reads 超过可用数量，按比例缩小
    if requested > total_reads:
        print(f"Warning: requested {requested:,} reads but only have {total_reads:,}. "
              f"Scaling down...", file=sys.stderr)
        sizes = (sizes * total_reads / requested).astype(int)
        sizes = np.maximum(sizes, 0)
        # 补齐因取整丢失的 reads
        diff = total_reads - int(sizes.sum())
        if diff > 0:
            active = np.where(sizes > 0)[0]
            if len(active) > 0:
                add_idx = np.random.choice(active, size=min(diff, len(active)), replace=False)
                sizes[add_idx] += 1

    # 分配并输出（有放回抽样）
    print(f"Writing {args.output}...", file=sys.stderr)
    opener = gzip.open if args.output.endswith('.gz') else open
    written = 0
    with opener(args.output, 'wt') as out:
        for b in range(args.n_barcodes):
            if sizes[b] == 0:
                continue
            barcode = f"{args.barcode_prefix}{b:04d}"
            for _ in range(sizes[b]):
                name, seq, plus, qual = random.choice(records)
                # 解析原始 read name，提取 {name}/{pair} 部分（去掉原有 barcode 和 UMI）
                # 原始格式: @ILLUMINA_NAME/1_barcode_UMI
                base_name = name.split('_')[0]
                # UMI 尾部加随机数，避免下游按 (barcode, UMI) 去重时被误合并
                umi_rand = ''.join(random.choices('ACGTN', k=6))
                new_name = f"{base_name}_{barcode}_{umi_rand}#empty"
                out.write(f"{new_name}\n{seq}\n{plus}\n{qual}\n")
                written += 1

    used_barcodes = np.sum(sizes > 0)
    print(f"  Written: {written:,} reads across {used_barcodes} barcodes", file=sys.stderr)
    print(f"  Barcode size: mean={sizes.mean():.1f}, std={sizes.std():.1f}, "
          f"median={np.median(sizes):.0f}, min={sizes.min()}, max={sizes.max()}",
          file=sys.stderr)


if __name__ == '__main__':
    main()
