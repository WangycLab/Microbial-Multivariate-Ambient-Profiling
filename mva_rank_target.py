#!/usr/bin/env python3
"""
从 counts.tsv 中按 UMI 排名选取目标 UMI 数：
  1) 按 UMI 降序排名（rank 1 = UMI 最高的 barcode）；
  2) 从 --min-rank 开始，合并所有 -m 值作为总抽样数，连续截取对应 barcode；
  3) 随机打乱后按各 -m 值依次分配到各 group，输出 UMI 数（每行一个数值，无表头）。
列约定：第 1 列 barcode，第 3 列 umi_count（与原脚本一致）。
"""

import argparse
import random
import sys


def main():
    parser = argparse.ArgumentParser(
        description="按 UMI 排名选取目标 UMI 数，-m 可多次传入以指定各 group 大小"
    )
    parser.add_argument("-i", "--input", required=True,
                        help="输入 counts.tsv（第1列 barcode，第3列 umi_count）")
    parser.add_argument("--min-rank", type=int, required=True,
                        help="起始排名，1-based，按 UMI 降序（rank 1 = UMI 最高）")
    parser.add_argument("-m", type=int, action="append", required=True,
                        help="每个 group 的 barcode 数，可多次传入（如 -m 500 -m 300）")
    parser.add_argument("-o", "--output", default="umi_targets.txt",
                        help="单 group 时的输出文件名（仅传入一个 -m 时生效）")
    parser.add_argument("--prefix", default="group",
                        help="多 group 输出前缀（传入多个 -m 时生效，如 group1.txt）")
    parser.add_argument("--seed", type=int, default=None,
                        help="随机种子，控制随机分配（多 group 时生效）")
    args = parser.parse_args()

    if args.min_rank < 1:
        print("[error] --min-rank 必须 >= 1", file=sys.stderr)
        sys.exit(1)
    for idx, m in enumerate(args.m):
        if m < 1:
            print(f"[error] 第 {idx + 1} 个 -m ({m}) 必须 >= 1", file=sys.stderr)
            sys.exit(1)

    # 读取 counts.tsv：第1列 barcode，第3列 umi_count
    records = []
    with open(args.input, "r") as f:
        header = f.readline()  # 跳过表头
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            bc = fields[0]
            try:
                umi = int(fields[2])
            except ValueError:
                continue
            records.append((bc, umi))

    if not records:
        print("[error] no barcodes found in input", file=sys.stderr)
        sys.exit(1)

    # 按 UMI 降序排名（rank 1 = UMI 最高的 barcode）
    records.sort(key=lambda r: r[1], reverse=True)
    n_total = len(records)
    print(f"[info] {n_total} barcodes loaded, sorted by UMI descending", file=sys.stderr)

    # 合并所有 -m 作为总抽样数，从 min-rank 起连续截取
    total_m = sum(args.m)
    start = args.min_rank - 1
    if start >= n_total:
        print(f"[error] --min-rank {args.min_rank} 超出 barcode 总数 {n_total}",
              file=sys.stderr)
        sys.exit(1)

    end = min(start + total_m, n_total)
    selected = [umi for _, umi in records[start:end]]
    actual_total = len(selected)
    if actual_total < total_m:
        print(f"[warn] 从 rank {args.min_rank} 起仅剩 {actual_total} 个 barcode"
              f"（请求 {total_m}），全部选取", file=sys.stderr)

    n_groups = len(args.m)

    if n_groups == 1:
        # 单 group：输出到 -o
        with open(args.output, "w") as f:
            for t in selected:
                f.write(f"{t}\n")
        print(f"[info] wrote {actual_total} UMI targets to {args.output}", file=sys.stderr)
    else:
        # 多 group：随机打乱后按各 -m 值依次分配
        rng = random.Random(args.seed)
        shuffled = selected[:]
        rng.shuffle(shuffled)

        idx = 0
        for g, m in enumerate(args.m):
            # 按比例缩放：实际可分配数不足时按 -m 比例分配
            if actual_total < total_m:
                size = round(m * actual_total / total_m) if g < n_groups - 1 \
                    else actual_total - idx
            else:
                size = m
            size = min(size, actual_total - idx)
            group_targets = shuffled[idx:idx + size]
            idx += size
            out_file = f"{args.prefix}{g + 1}.txt"
            with open(out_file, "w") as f:
                for t in group_targets:
                    f.write(f"{t}\n")
            print(f"[info] wrote {len(group_targets)} UMI targets to {out_file}",
                  file=sys.stderr)


if __name__ == "__main__":
    main()
