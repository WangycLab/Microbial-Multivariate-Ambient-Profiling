#!/usr/bin/env python3
"""
统计simulated.fq中每个barcode的reads数和UMI数（数据已去重）
Usage: python count_barcode_umi.py -i simulated.fq -o output.tsv
"""

import sys
import argparse
from collections import defaultdict

def parse_fastq_readname(readname):
    """
    从read name中提取barcode和UMI
    格式: @LH29866:259:22NVVYLT4:6:1115:3027:4879/1_barcode_UMI
    """
    # 去掉开头的@
    if readname.startswith('@'):
        readname = readname[1:]
    
    # 按照_分割
    parts = readname.split('_')
    
    if len(parts) >= 3:
        # 格式: .../1_barcode_UMI
        barcode = parts[-2]  # ACTTCCGACCCGATTGTATC#large
        umi = parts[-1]      # TACGTACC
        return barcode, umi
    else:
        return None, None

def count_barcode_umi(fastq_file):
    """
    统计每个barcode的reads数和UMI数
    数据已去重，所以umi_count = reads_count
    """
    barcode_reads = defaultdict(int)
    
    with open(fastq_file, 'r') as f:
        line_num = 0
        for line in f:
            line = line.strip()
            
            # FASTQ格式: 每4行一个read
            if line_num % 4 == 0:
                # read name行
                if line.startswith('@'):
                    barcode, umi = parse_fastq_readname(line)
                    if barcode:
                        barcode_reads[barcode] += 1
            line_num += 1
    
    return barcode_reads

def main():
    parser = argparse.ArgumentParser(description='统计barcode的reads数和UMI数（数据已去重）')
    parser.add_argument('-i', '--input', required=True, 
                        help='输入的FASTQ文件路径')
    parser.add_argument('-o', '--output', required=True,
                        help='输出的TSV文件路径')
    
    args = parser.parse_args()
    
    print(f"正在处理: {args.input}", file=sys.stderr)
    
    barcode_reads = count_barcode_umi(args.input)
    
    # 按 umi_counts 降序排序
    sorted_barcodes = sorted(barcode_reads.items(), key=lambda x: x[1], reverse=True)
    
    # 写入输出文件
    with open(args.output, 'w') as f:
        f.write("barcode\treads_count\tumi_count\n")
        for barcode, count in sorted_barcodes:
            # 数据已去重，umi_count = reads_count
            f.write(f"{barcode}\t{count}\t{count}\n")
    
    print(f"完成！共统计 {len(barcode_reads)} 个barcode", file=sys.stderr)
    print(f"输出文件: {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()