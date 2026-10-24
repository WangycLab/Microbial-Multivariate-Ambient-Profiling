#!/usr/bin/env python3
"""
UMI dedup: one-pass BAM → dedup.bam + out.fq + counts.tsv

第一遍读 BAM，遇到新的 (contig, gene, UMI) key 时：
  - 直接写 dedup.bam
  - 直接从 BAM record 提取序列/质量写 out.fq
无需 keep_reads.txt 中间过滤，无需 seqkit grep。
"""

import argparse
import pysam
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-b", "--bam", required=True, help="sorted BAM file")
parser.add_argument("-c", "--counts", required=True, help="counts.tsv")
parser.add_argument("-d", "--dedup-bam", required=True, help="output dedup BAM")
parser.add_argument("-o", "--outfq", required=True, help="output FASTQ")
parser.add_argument("-k", "--keep", default="/dev/null", help="keep reads list (for record)")
parser.add_argument("-g", "--gene-delimiter", default="#", help="RNAME 中 contig 与 gene 的分隔符; 默认 #")
args = parser.parse_args()

read_count = defaultdict(int)
umi_sets = defaultdict(set)
keep_reads = set()

bam_in = pysam.AlignmentFile(args.bam, "rb")
bam_out = pysam.AlignmentFile(args.dedup_bam, "wb", header=bam_in.header)
fq_out = open(args.outfq, "w")

for r in bam_in:
    if r.is_unmapped or r.is_secondary or r.is_supplementary:
        continue

    qname = r.query_name
    try:
        barcode = qname.split("_")[1]
        umi = qname.split("_")[2]
    except Exception:
        continue

    rname = bam_in.get_reference_name(r.reference_id)
    if args.gene_delimiter in rname:
        ref_genome, gene = rname.split(args.gene_delimiter, 1)
    else:
        ref_genome, gene = rname, rname

    read_count[barcode] += 1
    key = (ref_genome, gene, umi)

    if key not in umi_sets[barcode]:
        umi_sets[barcode].add(key)
        keep_reads.add(qname)

        # 直接写 dedup.bam
        bam_out.write(r)

        # 从 BAM record 提取序列/质量，直接写 FASTQ
        seq = r.query_sequence
        qual_arr = r.query_qualities
        if seq and qual_arr and len(seq) == len(qual_arr):
            qual_str = "".join(chr(q + 33) for q in qual_arr)
            fq_out.write(f"@{qname}\n{seq}\n+\n{qual_str}\n")

bam_in.close()
bam_out.close()
fq_out.close()

# counts.tsv
with open(args.counts, "w") as out:
    out.write("barcode\treads_count\tumi_count\n")
    for bc in sorted(read_count, key=lambda x: len(umi_sets[x]), reverse=True):
        out.write(f"{bc}\t{read_count[bc]}\t{len(umi_sets[bc])}\n")

# keep_reads.txt (可选，仅作记录)
with open(args.keep, "w") as f:
    for r in keep_reads:
        f.write(r + "\n")

print(f"[INFO] total barcodes: {len(read_count)}")
print(f"[INFO] unique kept reads: {len(keep_reads)}")
print(f"[INFO] dedup BAM → {args.dedup_bam}")
print(f"[INFO] dedup FASTQ → {args.outfq}")
print("[DONE]")
