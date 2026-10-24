#!/usr/bin/env python3
"""
Cell Ranger 风格 knee point 检测脚本（精简输出版）。

算法：Kneedle / 最大曲率法
1. 按 UMI count 降序排列 barcode
2. 计算归一化累积 UMI 曲线
3. 找到曲线到首尾连线距离最大的点，即为 knee point

输出：仅一个数字（默认 n_cells，即过滤后保留的 barcode 数）。
"""

import sys
import argparse
import numpy as np


def find_knee_point(umi_counts: np.ndarray, sensitivity: float = 1.0, smooth: int = 1) -> dict:
    """
    对给定的 UMI count 数组（降序），用 Kneedle 算法找到 knee point。

    sensitivity: Kneedle 敏感度系数（S），控制检测严格程度。
                 值越大 knee 越靠左（更严格），值越小越宽松。
    smooth:      平滑窗口大小，对归一化累积曲线做移动平均以减少噪声。
                 1 表示不平滑。

    返回包含 knee index、UMI 阈值、过滤后 barcode 数量等信息的字典。
    """
    n = len(umi_counts)
    if n < 3:
        return {"knee_index": 0, "umi_threshold": 0, "n_cells": n, "n_total": n}

    # 累积 UMI
    cumulative = np.cumsum(umi_counts, dtype=np.float64)
    total_umi = cumulative[-1]

    # 归一化到 [0, 1]
    x = np.arange(1, n + 1) / n          # 归一化 rank
    y = cumulative / total_umi            # 归一化累积 UMI

    # 平滑处理
    if smooth > 1 and n >= smooth:
        kernel = np.ones(smooth) / smooth
        y = np.convolve(y, kernel, mode="same")

    # 首尾连线: y = x（在归一化空间中）
    # 点到线的距离 = |y - x| / sqrt(2)
    distances = np.abs(y - x) / np.sqrt(2)

    # Kneedle 算法: 从右到左扫描，找到满足 sensitivity 条件的 knee
    d_max = np.max(distances)
    d_min = np.min(distances)
    threshold = d_min + sensitivity * (d_max - d_min)

    candidates = np.where(distances > threshold)[0]
    if len(candidates) == 0:
        knee_idx = int(np.argmax(distances))
    else:
        knee_idx = int(candidates[0])

    result = {
        "knee_index": knee_idx,
        "umi_threshold": int(umi_counts[knee_idx]),
        "n_cells": knee_idx + 1,                # 过滤后保留的 barcode 数
        "n_total": n,
        "cumulative_umi_at_knee": float(cumulative[knee_idx]),
        "total_umi": float(total_umi),
        "fraction_umi": float(cumulative[knee_idx] / total_umi),
        "fraction_barcodes": float((knee_idx + 1) / n),
    }
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Cell Ranger 风格 knee point 检测（输出单个数字）"
    )
    parser.add_argument(
        "input",
        nargs="?",
        default="simulated_counts.tsv",
        help="输入 TSV 文件路径 (默认: simulated_counts.tsv)",
    )
    parser.add_argument(
        "--library",
        choices=["large", "small", "all"],
        default="all",
        help="选择处理的文库: large, small, 或 all (默认: all)",
    )
    parser.add_argument(
        "--column",
        default="umi_count",
        help="用于排序的列名 (默认: umi_count)",
    )
    parser.add_argument(
        "--sensitivity",
        type=float,
        default=1.0,
        help="Kneedle 敏感度系数 S (默认: 1.0)。值越大 knee 越靠左(更严格)，越小越宽松",
    )
    parser.add_argument(
        "--smooth",
        type=int,
        default=1,
        help="平滑窗口大小 (默认: 1，不平滑)。对归一化累积曲线做移动平均以减少噪声",
    )
    parser.add_argument(
        "--value",
        choices=["n_cells", "umi_threshold", "knee_index", "fraction_umi", "fraction_barcodes"],
        default="n_cells",
        help="输出哪个值: n_cells(细胞数), umi_threshold(UMI阈值), "
             "knee_index(rank), fraction_umi(UMI占比), fraction_barcodes(barcode占比) "
             "(默认: n_cells)",
    )
    args = parser.parse_args()

    # 读取 TSV
    try:
        with open(args.input) as f:
            header = f.readline().strip().split("\t")
            lines = f.readlines()
    except FileNotFoundError:
        print(f"错误: 找不到文件 '{args.input}'", file=sys.stderr)
        sys.exit(1)

    # 解析数据
    col_idx = header.index(args.column)

    barcodes = []
    counts = []
    for line in lines:
        if not line.strip():
            continue
        fields = line.strip().split("\t")
        barcode = fields[0]
        umi = int(fields[col_idx])

        # 按文库过滤
        if args.library == "large" and "#large" not in barcode:
            continue
        if args.library == "small" and "#small" not in barcode:
            continue

        barcodes.append(barcode)
        counts.append(umi)

    if not counts:
        print("错误: 没有匹配的数据", file=sys.stderr)
        sys.exit(1)

    # 按 UMI count 降序排列
    order = np.argsort(counts)[::-1]
    umi_sorted = np.array([counts[i] for i in order], dtype=np.int64)

    # 检测 knee point
    result = find_knee_point(umi_sorted, sensitivity=args.sensitivity, smooth=args.smooth)

    # 只输出一个数字
    val = result[args.value]
    if args.value in ("fraction_umi", "fraction_barcodes"):
        print(f"{val:.6f}")
    else:
        print(int(val))


if __name__ == "__main__":
    main()
