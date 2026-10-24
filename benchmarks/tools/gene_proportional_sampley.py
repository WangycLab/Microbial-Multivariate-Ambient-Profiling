#!/usr/bin/env python3
"""
gene_proportional_sample.py
===========================
从 gene-annotated BAM 中按等比例抽取 reads，输出到 FASTQ。

分类抽样策略:
  - rRNA gene: 直接按 gene-id 抽样
  - non-rRNA gene: 先映射到 KEGG Module，按 module 抽样
    （一个 gene 对应多个 module 时，read 同时计入所有 module）
    无 module 映射的 non-rRNA reads 不参与抽样

用法:
  python gene_proportional_sample.py -i dedup.fq -b gene.sorted.bam \
      -g annotation.gff --gene-module gene_module.tsv -o sampled.fq \
      --coef 0.5 --n-small 2000

流程:
  1. 加载 rRNA gene set + gene->module 映射
  2. 扫描 BAM，按 barcode 分类:
     - rRNA reads → bc_gene_reads[bc][gene]
     - non-rRNA mapped reads → bc_mod_reads[bc][module]
     - non-rRNA unmapped → 跳过
  3. 按 total read count 降序取前 n_small 个 barcode
  4. rRNA 按 gene 概率比例抽样；non-rRNA 按 module 概率比例抽样
  5. 扫描 FASTQ，输出抽中的 reads，header barcode 后加 #small
"""

import sys
import gzip
import math
import argparse
import subprocess
import random
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description='Proportional sampling: rRNA by gene, non-rRNA by KEGG module')
    p.add_argument('-i', '--fq', required=True, help='Input FASTQ')
    p.add_argument('-b', '--bam', required=True, help='Gene-annotated sorted BAM')
    p.add_argument('-g', '--gff', required=True, help='GFF3 for rRNA gene identification')
    p.add_argument('--gene-module', required=True, help='gene to KEGG Module TSV (emapper output)')
    p.add_argument('-o', '--output', default='sampled.fq', help='Output FASTQ')
    p.add_argument('--coef-min', type=float, default=0.2,
                   help='Sampling coefficient lower bound. Default: 0.2')
    p.add_argument('--coef-max', type=float, default=0.5,
                   help='Sampling coefficient upper bound. Default: 0.5')
    p.add_argument('--n-small', type=int, default=2000,
                   help='Number of barcodes to sample (with replacement). Default: 2000')
    p.add_argument('--seed', type=int, default=42, help='Random seed')
    p.add_argument('--threads', type=int, default=8, help='samtools threads')
    return p.parse_args()


def load_rrna_genes(gff_path: str) -> set:
    rrna = set()
    print(f'[gff] loading rRNA genes from {gff_path}', file=sys.stderr)
    with open(gff_path, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            attrs = {}
            for raw_attr in fields[8].split(';'):
                attr = raw_attr.strip()
                if not attr or '=' not in attr:
                    continue
                key, value = attr.split('=', 1)
                attrs[key] = value
            if attrs.get('gene_biotype') == 'rRNA' and attrs.get('ID'):
                rrna.add(attrs['ID'])
    print(f'[gff] rRNA genes={len(rrna):,}', file=sys.stderr)
    return rrna


def load_gene_to_module(path: str) -> dict:
    """加载 gene->module 映射。格式: #query\tKEGG_Module"""
    gene_to_module = {}
    print(f'[module] loading gene->module map from {path}', file=sys.stderr)
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            gene, mod_text = parts[0], parts[1]
            if mod_text == '-' or not mod_text:
                continue
            modules = [m.strip() for m in mod_text.split(',') if m.strip() and m.strip() != '-']
            if modules:
                gene_to_module[gene] = modules
    n_assignments = sum(len(v) for v in gene_to_module.values())
    print(f'[module] mapped genes={len(gene_to_module):,} module assignments={n_assignments:,}',
          file=sys.stderr)
    return gene_to_module


def parse_bam(bam_path, rrna_genes, gene_to_module, threads=8):
    """
    扫描 BAM，分类返回:
      bc_rrna_reads: {barcode: {gene: [read_names]}}  — rRNA 按 gene-id
      bc_mod_reads:  {barcode: {module: [read_names]}} — non-rRNA 按 module
      bc_total:      {barcode: total_reads}
    """
    cmd = ['samtools', 'view', '-@', str(threads), bam_path]
    bc_rrna_reads = defaultdict(lambda: defaultdict(list))
    bc_mod_reads = defaultdict(lambda: defaultdict(list))
    bc_total = defaultdict(int)

    total = 0
    mapped = 0
    n_rrna = 0
    n_nonrrna_mapped = 0
    n_nonrrna_no_module = 0

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    for line in proc.stdout:
        total += 1
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 10:
            continue
        if int(fields[1]) & 4:  # unmapped
            continue

        ref = fields[2]
        if '#' in ref:
            gene = ref.split('#', 1)[1]
        else:
            gene = ref

        qname = fields[0]

        parts = qname.split('_')
        if len(parts) < 3:
            continue
        bc = parts[-2]

        bc_total[bc] += 1
        mapped += 1

        if gene in rrna_genes:
            # rRNA: 按 gene-id 归类
            bc_rrna_reads[bc][gene].append(qname)
            n_rrna += 1
        else:
            # non-rRNA: 映射到 module
            modules = gene_to_module.get(gene)
            if not modules:
                # 无 module 映射：按 gene-id 归类，与 rRNA 一起按 gene-id 抽样
                bc_rrna_reads[bc][gene].append(qname)
                n_nonrrna_no_module += 1
                continue
            for mod in modules:
                bc_mod_reads[bc][mod].append(qname)
            n_nonrrna_mapped += 1

        if total % 5000000 == 0:
            print(f'  ... BAM records: {total:,}', file=sys.stderr)

    proc.wait()
    print(f'  BAM total: {total:,}, Mapped: {mapped:,}', file=sys.stderr)
    print(f'    rRNA: {n_rrna:,}, non-rRNA with module: {n_nonrrna_mapped:,}, '
          f'non-rRNA no module: {n_nonrrna_no_module:,}', file=sys.stderr)
    return bc_rrna_reads, bc_mod_reads, bc_total


def probabilistic_sample(reads, coef):
    """概率比例抽样：floor + 小数部分按概率多抽 1 个。"""
    expected = len(reads) * coef
    n_keep = math.floor(expected)
    if random.random() < (expected - n_keep):
        n_keep += 1
    return min(n_keep, len(reads))


def main():
    args = parse_args()
    random.seed(args.seed)

    # === 1. 加载 rRNA genes + gene-to-module 映射 ===
    rrna_genes = load_rrna_genes(args.gff)
    gene_to_module = load_gene_to_module(args.gene_module)

    # === 2. 扫描 BAM ===
    print('[1/3] Parsing BAM...', file=sys.stderr)
    bc_rrna_reads, bc_mod_reads, bc_total = parse_bam(
        args.bam, rrna_genes, gene_to_module, threads=args.threads
    )

    # === 3. 有放回抽样 n_small 个 barcode ===
    print('\n[2/3] Selecting barcodes (with replacement)...', file=sys.stderr)
    all_bcs = list(bc_total.keys())
    umi_min = min(bc_total.values())
    umi_max = max(bc_total.values())
    umi_range = umi_max - umi_min if umi_max > umi_min else 1.0

    def umi_to_coef(umi):
        """UMI 越大 -> coef 越小 (downsampling 更狠); UMI 越小 -> coef 越大"""
        frac = (umi - umi_min) / umi_range
        return args.coef_max - frac * (args.coef_max - args.coef_min)

    selected = []  # [(orig_bc, new_bc, coef), ...]
    for _ in range(args.n_small):
        orig_bc = random.choice(all_bcs)
        new_bc = ''.join(random.choices('ATCG', k=10))
        coef = umi_to_coef(bc_total[orig_bc])
        selected.append((orig_bc, new_bc, coef))
    print(f'  Selected {len(selected):,} barcodes (with replacement)', file=sys.stderr)

    # === 4. 分层抽样 ===
    print(f'\n[3/3] Proportional sampling (coef range {args.coef_min}-{args.coef_max})...', file=sys.stderr)
    read_to_new_bcs = defaultdict(list)  # qname -> [new_bc1, new_bc2, ...]
    n_rrna_keep = 0
    n_mod_keep = 0

    for orig_bc, new_bc, coef in selected:
        # rRNA: 按 gene-id 抽样
        for gene, reads in bc_rrna_reads[orig_bc].items():
            n_keep = probabilistic_sample(reads, coef)
            if n_keep > 0:
                for r in random.sample(reads, n_keep):
                    read_to_new_bcs[r].append(new_bc)
                    n_rrna_keep += 1

        # non-rRNA: 按 module 抽样
        for mod, reads in bc_mod_reads[orig_bc].items():
            n_keep = probabilistic_sample(reads, coef)
            if n_keep > 0:
                for r in random.sample(reads, n_keep):
                    read_to_new_bcs[r].append(new_bc)
                    n_mod_keep += 1

    print(f'  rRNA reads to keep: {n_rrna_keep:,}', file=sys.stderr)
    print(f'  non-rRNA reads to keep: {n_mod_keep:,}', file=sys.stderr)
    print(f'  Unique read names: {len(read_to_new_bcs):,}', file=sys.stderr)

    # === 5. 扫描 FASTQ 输出 ===
    print('\n[4/4] Writing output FASTQ...', file=sys.stderr)
    opener = gzip.open if args.fq.endswith('.gz') else open
    kept = 0

    with opener(args.fq, 'rt') as fin, open(args.output, 'w') as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()

            if not header.startswith('@'):
                continue

            qname = header.rstrip('\n')[1:]  # 去掉 @
            if qname in read_to_new_bcs:
                for new_bc in read_to_new_bcs[qname]:
                    parts = qname.split('_')
                    if len(parts) >= 3:
                        parts[-2] = f"{new_bc}#small"
                        new_header = '@' + '_'.join(parts) + '\n'
                    else:
                        new_header = header
                    fout.write(new_header)
                    fout.write(seq)
                    fout.write(plus)
                    fout.write(qual)
                    kept += 1

    print(f'  Output reads: {kept:,}', file=sys.stderr)
    print(f'\nDone! {args.output}', file=sys.stderr)


if __name__ == '__main__':
    main()
