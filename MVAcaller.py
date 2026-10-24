#!/usr/bin/env python3
"""
MVAcaller (multivariate ambient caller): cell calling for microbial single-cell RNA-seq.

Implements a three-test strategy that combines:
  1. rRNA-fraction test   — barcodes with rRNA proportion inconsistent with ambient
  2. KEGG-module test     — barcodes whose KO composition deviates from ambient
  3. non-rRNA (or total-RNA) composition test  — EmptyDrops-style multinomial test on the
     protein-coding gene profile (or total transcriptome in --mode-omics)

P-values from the three tests are combined via Fisher's method, followed by
Benjamini-Hochberg FDR correction. Barcodes with adjusted combined P < --pvalue-thre
(default 0.01) are called as cells.

Two input modes are supported:
  * Memory mode:    parse BAM/SAM directly into memory (simpler, smaller datasets)
  * SQLite mode:    build a SQLite molecule DB first, then query (fast for large BAMs,
                    reusable across runs)

Quick start (memory mode):
    python MVAcaller.py \\
        -i gene-tagged.bam \\
        -g annotation.gff3 \\
        --kegg-map gene_to_kegg.tsv \\
        -o output/prefix \\
        --qname-simple \\
        --umi-lower 50 --umi-upper 5000

Quick start (SQLite mode, reusable):
    python MVAcaller.py \\
        -i gene-tagged.bam \\
        -g annotation.gff3 \\
        --kegg-map gene_to_kegg.tsv \\
        --build-sqlite molecules.sqlite --reuse-sqlite \\
        -o output/prefix \\
        --qname-simple
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import math
import multiprocessing as mp
import os
import random
import re
import sqlite3
import sys
import time
from collections import Counter, defaultdict
from typing import Dict, Iterator, List, Optional, Sequence, Tuple


CellFeatureCounts = Dict[str, Counter]


# ---------------------------------------------------------------------------
# BAM / SAM 解析
# ---------------------------------------------------------------------------

def smart_open_text(path: str):
    if path == "-":
        return sys.stdin
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def iter_sam_fields(path: str) -> Iterator[List[str]]:
    is_probably_bam = path != "-" and path.lower().endswith((".bam", ".cram"))
    if is_probably_bam:
        try:
            import pysam  # type: ignore
        except ImportError:
            raise SystemExit("BAM/CRAM input requires pysam: pip install pysam") from None
        mode = "rc" if path.lower().endswith(".cram") else "rb"
        with pysam.AlignmentFile(path, mode, ignore_truncation=True) as bam:
            for read in bam:
                yield read.to_string().split("\t")
        return

    with smart_open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 11:
                yield fields


def parse_cb_umi_from_qname(qname: str, regex: re.Pattern[str]) -> Optional[Tuple[str, str]]:
    match = regex.search(qname)
    if not match:
        return None
    return match.group("cb"), match.group("umi")


def parse_cb_umi_from_tags(fields: List[str]) -> Optional[Tuple[str, str]]:
    cb = None
    umi = None
    for field in fields[11:]:
        if field.startswith("CB:Z:"):
            cb = field[5:]
        elif field.startswith("UB:Z:"):
            umi = field[5:]
        elif field.startswith("UR:Z:"):
            umi = field[5:]
    if cb is None or umi is None:
        return None
    return cb, umi


def gene_from_rname(rname: str, delimiter: str) -> Optional[str]:
    if rname == "*" or delimiter not in rname:
        return None
    gene = rname.rsplit(delimiter, 1)[1]
    return gene or None


# ---------------------------------------------------------------------------
# 内存模式：直接从 BAM 加载到内存
# ---------------------------------------------------------------------------

def build_counts_from_bam(
    bam_path: str,
    qname_regex: Optional[re.Pattern[str]],
    gene_delimiter: str,
    min_mapq: int,
    skip_secondary: bool,
    progress_every: int,
    use_cb_tag: bool = False,
) -> Tuple[CellFeatureCounts, Counter, Counter]:
    """Parse BAM and return (cell_gene_counts, umi_per_cb, stats). No deduplication."""
    stats = Counter()
    started = time.time()

    cell_gene_counts: CellFeatureCounts = defaultdict(Counter)
    umi_per_cb = Counter()
    genes_seen = set()

    for fields in iter_sam_fields(bam_path):
        stats["reads_total"] += 1
        qname = fields[0]
        try:
            flag = int(fields[1])
            mapq = int(fields[4])
        except ValueError:
            stats["reads_bad_sam"] += 1
            continue

        if flag & 0x4:
            stats["reads_unmapped"] += 1
            continue
        if skip_secondary and (flag & 0x100 or flag & 0x800):
            stats["reads_secondary_or_supplementary"] += 1
            continue
        if mapq < min_mapq:
            stats["reads_low_mapq"] += 1
            continue

        gene = gene_from_rname(fields[2], gene_delimiter)
        if gene is None:
            stats["reads_no_gene"] += 1
            continue

        if use_cb_tag:
            cb_umi = parse_cb_umi_from_tags(fields)
        else:
            cb_umi = parse_cb_umi_from_qname(qname, qname_regex) if qname_regex is not None else None
        if cb_umi is None:
            stats["reads_no_cb_umi"] += 1
            continue

        cb, umi = cb_umi

        cell_gene_counts[cb][gene] += 1
        umi_per_cb[cb] += 1
        genes_seen.add(gene)
        stats["reads_counted"] += 1

        if progress_every and stats["reads_total"] % progress_every == 0:
            elapsed = max(time.time() - started, 1e-6)
            print(
                f"[bam] reads={stats['reads_total']:,} counted={stats['reads_counted']:,} "
                f"rate={stats['reads_total']/elapsed:,.0f}/s",
                file=sys.stderr,
                flush=True,
            )

    stats["molecules_total"] = stats["reads_counted"]
    stats["barcodes_detected"] = len(umi_per_cb)
    stats["genes_detected"] = len(genes_seen)

    elapsed = max(time.time() - started, 1e-6)
    print(
        f"[bam] done reads={stats['reads_total']:,} counted={stats['reads_counted']:,} "
        f"molecules={stats['molecules_total']:,} barcodes={stats['barcodes_detected']:,} "
        f"genes={stats['genes_detected']:,} time={elapsed:.1f}s",
        file=sys.stderr,
        flush=True,
    )
    return cell_gene_counts, umi_per_cb, stats


# ---------------------------------------------------------------------------
# SQLite 模式
# ---------------------------------------------------------------------------

def build_sqlite(
    bam_path: str,
    qname_regex: Optional[re.Pattern[str]],
    gene_delimiter: str,
    min_mapq: int,
    skip_secondary: bool,
    sqlite_path: str,
    batch_size: int,
    progress_every: int,
    use_cb_tag: bool = False,
) -> Tuple[CellFeatureCounts, Counter, Counter]:
    os.makedirs(os.path.dirname(os.path.abspath(sqlite_path)) or ".", exist_ok=True)

    conn = sqlite3.connect(sqlite_path)
    cur = conn.cursor()
    cur.execute("PRAGMA journal_mode=WAL")
    cur.execute("PRAGMA synchronous=OFF")
    cur.execute("PRAGMA temp_store=MEMORY")
    cur.execute("PRAGMA cache_size=-200000")
    cur.execute("DROP TABLE IF EXISTS molecules")
    cur.execute(
        "CREATE TABLE molecules (cb TEXT NOT NULL, umi TEXT NOT NULL, gene TEXT NOT NULL)"
    )
    cur.execute("CREATE INDEX idx_molecules_cb ON molecules(cb)")
    cur.execute("CREATE INDEX idx_molecules_cb_gene ON molecules(cb, gene)")
    conn.commit()

    stats = Counter()
    batch = []
    started = time.time()

    for fields in iter_sam_fields(bam_path):
        stats["reads_total"] += 1
        qname = fields[0]
        try:
            flag = int(fields[1])
            mapq = int(fields[4])
        except ValueError:
            stats["reads_bad_sam"] += 1
            continue

        if flag & 0x4:
            stats["reads_unmapped"] += 1
            continue
        if skip_secondary and (flag & 0x100 or flag & 0x800):
            stats["reads_secondary_or_supplementary"] += 1
            continue
        if mapq < min_mapq:
            stats["reads_low_mapq"] += 1
            continue

        gene = gene_from_rname(fields[2], gene_delimiter)
        if gene is None:
            stats["reads_no_gene"] += 1
            continue

        if use_cb_tag:
            cb_umi = parse_cb_umi_from_tags(fields)
        else:
            cb_umi = parse_cb_umi_from_qname(qname, qname_regex) if qname_regex is not None else None
        if cb_umi is None:
            stats["reads_no_cb_umi"] += 1
            continue

        cb, umi = cb_umi

        batch.append((cb, umi, gene))
        stats["reads_counted"] += 1

        if len(batch) >= batch_size:
            cur.executemany("INSERT INTO molecules VALUES (?, ?, ?)", batch)
            conn.commit()
            batch.clear()

        if progress_every and stats["reads_total"] % progress_every == 0:
            elapsed = max(time.time() - started, 1e-6)
            print(
                f"[sqlite] reads={stats['reads_total']:,} counted={stats['reads_counted']:,} "
                f"rate={stats['reads_total']/elapsed:,.0f}/s",
                file=sys.stderr,
                flush=True,
            )

    if batch:
        cur.executemany("INSERT INTO molecules VALUES (?, ?, ?)", batch)
        conn.commit()

    cell_gene_counts: CellFeatureCounts = defaultdict(Counter)
    genes_seen = set()
    for cb, gene, count in cur.execute("SELECT cb, gene, COUNT(*) FROM molecules GROUP BY cb, gene"):
        cell_gene_counts[cb][gene] = count
        genes_seen.add(gene)

    umi_per_cb = Counter()
    for cb, count in cur.execute("SELECT cb, COUNT(*) FROM molecules GROUP BY cb"):
        umi_per_cb[cb] = count

    row = cur.execute("SELECT COUNT(*) FROM molecules").fetchone()
    stats["molecules_total"] = int(row[0])
    stats["barcodes_detected"] = len(umi_per_cb)
    stats["genes_detected"] = len(genes_seen)
    conn.close()

    elapsed = max(time.time() - started, 1e-6)
    print(
        f"[sqlite] done reads={stats['reads_total']:,} counted={stats['reads_counted']:,} "
        f"molecules={stats['molecules_total']:,} barcodes={stats['barcodes_detected']:,} "
        f"genes={stats['genes_detected']:,} time={elapsed:.1f}s",
        file=sys.stderr,
        flush=True,
    )
    return cell_gene_counts, umi_per_cb, stats


def load_gene_counts_sqlite(sqlite_path: str) -> Tuple[CellFeatureCounts, Counter, Counter]:
    if not os.path.exists(sqlite_path):
        raise SystemExit(f"SQLite molecule DB does not exist: {sqlite_path}")
    print(f"[count:sqlite] loading {sqlite_path}", file=sys.stderr, flush=True)
    conn = sqlite3.connect(sqlite_path)
    cur = conn.cursor()

    cell_gene_counts: CellFeatureCounts = defaultdict(Counter)
    genes_seen = set()
    for cb, gene, count in cur.execute("SELECT cb, gene, COUNT(*) FROM molecules GROUP BY cb, gene"):
        cell_gene_counts[cb][gene] = count
        genes_seen.add(gene)

    umi_per_cb = Counter()
    for cb, count in cur.execute("SELECT cb, COUNT(*) FROM molecules GROUP BY cb"):
        umi_per_cb[cb] = count

    row = cur.execute("SELECT COUNT(*) FROM molecules").fetchone()
    conn.close()

    stats = Counter()
    stats["molecules_total"] = int(row[0])
    stats["barcodes_detected"] = len(umi_per_cb)
    stats["genes_detected"] = len(genes_seen)
    print(
        f"[count:sqlite] loaded molecules={stats['molecules_total']:,} "
        f"barcodes={stats['barcodes_detected']:,} genes={stats['genes_detected']:,}",
        file=sys.stderr,
        flush=True,
    )
    return cell_gene_counts, umi_per_cb, stats


# ---------------------------------------------------------------------------
# 注释加载
# ---------------------------------------------------------------------------

def parse_gff_attrs(attr_text: str) -> Dict[str, str]:
    attrs = {}
    for raw_attr in attr_text.split(";"):
        attr = raw_attr.strip()
        if not attr or "=" not in attr:
            continue
        key, value = attr.split("=", 1)
        attrs[key] = value
    return attrs


def load_rrna_genes(gff_path: str) -> set:
    rrna = set()
    print(f"[gff] loading rRNA genes from {gff_path}", file=sys.stderr, flush=True)
    with open(gff_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = parse_gff_attrs(fields[8])
            if attrs.get("gene_biotype") == "rRNA" and attrs.get("ID"):
                rrna.add(attrs["ID"])
    print(f"[gff] rRNA genes={len(rrna):,}", file=sys.stderr, flush=True)
    return rrna


def load_gene_to_ko(path: str, multi_ko: str) -> Dict[str, List[str]]:
    gene_to_ko: Dict[str, List[str]] = {}
    print(f"[kegg] loading gene->KO map from {path}", file=sys.stderr, flush=True)
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("##"):
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            gene, ko_text = fields[0], fields[1]
            if ko_text == "-" or not ko_text:
                continue
            kos = [ko.strip() for ko in ko_text.split(",") if ko.strip() and ko.strip() != "-"]
            if not kos:
                continue
            if multi_ko == "first":
                kos = kos[:1]
            elif multi_ko == "drop" and len(kos) > 1:
                continue
            gene_to_ko[gene] = kos
    n_assignments = sum(len(v) for v in gene_to_ko.values())
    print(
        f"[kegg] mapped genes={len(gene_to_ko):,} KO assignments={n_assignments:,} multi_ko={multi_ko}",
        file=sys.stderr,
        flush=True,
    )
    return gene_to_ko


# ---------------------------------------------------------------------------
# KO 计数转换
# ---------------------------------------------------------------------------

def gene_counts_to_ko_counts(cell_gene_counts: CellFeatureCounts, gene_to_ko: Dict[str, List[str]], rrna_genes: Optional[set] = None) -> Tuple[CellFeatureCounts, CellFeatureCounts, Counter]:
    """Split non-rRNA gene counts into KO-mapped and all non-rRNA gene counts per barcode.

    cell_ko_counts: non-rRNA genes mapped to KOs (split by KO).
    cell_non_rrna_counts: ALL non-rRNA genes (regardless of KO mapping), for non-rRNA EmptyDrops.
    """
    cell_ko_counts: CellFeatureCounts = defaultdict(Counter)
    cell_non_rrna_counts: CellFeatureCounts = defaultdict(Counter)
    stats = Counter()
    for cb, gene_counts in cell_gene_counts.items():
        for gene, count in gene_counts.items():
            if rrna_genes and gene in rrna_genes:
                stats["rrna_gene_umis_filtered"] += count
                continue
            # All non-rRNA genes go into cell_non_rrna_counts
            cell_non_rrna_counts[cb][gene] += count
            # Also split KO-mapped genes into cell_ko_counts
            kos = gene_to_ko.get(gene)
            if not kos:
                stats["non_rrna_no_ko_gene_umis"] += count
                continue
            for ko in kos:
                cell_ko_counts[cb][ko] += count
                stats["ko_feature_umis"] += count
            stats["mapped_gene_umis"] += count
    stats["barcodes_with_ko"] = len(cell_ko_counts)
    stats["ko_features_detected"] = len({ko for counts in cell_ko_counts.values() for ko in counts})
    stats["non_rrna_genes_detected"] = len({g for counts in cell_non_rrna_counts.values() for g in counts})
    return cell_ko_counts, cell_non_rrna_counts, stats


# ---------------------------------------------------------------------------
# 统计检验核心函数
# ---------------------------------------------------------------------------

def bh_adjust(pairs: List[Tuple[str, float]]) -> Dict[str, float]:
    if not pairs:
        return {}
    ordered = sorted(pairs, key=lambda item: item[1])
    n = len(ordered)
    adjusted = {}
    prev = 1.0
    for rank_from_end, (cb, pval) in enumerate(reversed(ordered), 1):
        rank = n - rank_from_end + 1
        qval = min(prev, pval * n / rank)
        adjusted[cb] = min(qval, 1.0)
        prev = qval
    return adjusted


def log_multinomial(counts: Counter, log_probs: Dict[str, float], log_factorial: List[float], floor_logp: float) -> float:
    total = sum(counts.values())
    value = log_factorial[total]
    for feature, count in counts.items():
        value -= log_factorial[count]
        value += count * log_probs.get(feature, floor_logp)
    return value


def make_log_factorial(n: int) -> List[float]:
    out = [0.0] * (n + 1)
    for i in range(2, n + 1):
        out[i] = out[i - 1] + math.log(i)
    return out


def simulate_lower_counts_chunk(
    args: Tuple[int, int, int, Sequence[float], Sequence[float], Sequence[int], Sequence[float]]
) -> List[int]:
    start_sim, n_sims, seed, probs, log_probs, depths, obs_by_candidate = args
    max_count = max(depths)
    cumulative = []
    total = 0.0
    for p in probs:
        total += p
        cumulative.append(total)

    depth_to_candidate_indices: Dict[int, List[int]] = defaultdict(list)
    for idx, depth in enumerate(depths):
        depth_to_candidate_indices[depth].append(idx)

    lower = [0] * len(depths)
    for isim in range(start_sim, start_sim + n_sims):
        rng = random.Random(seed * (isim + 1))
        curr_counts = [0] * len(probs)
        logp = 0.0
        for depth in range(1, max_count + 1):
            idx = bisect.bisect_left(cumulative, rng.random() * total)
            curr_counts[idx] += 1
            logp += log_probs[idx] + math.log(depth) - math.log(curr_counts[idx])
            for cand_idx in depth_to_candidate_indices.get(depth, ()):
                if logp < obs_by_candidate[cand_idx]:
                    lower[cand_idx] += 1
    return lower


def simulate_pvalues(
    candidates: List[str],
    depths: List[int],
    obs: Dict[str, float],
    sim_n: int,
    probs: List[float],
    log_probs: List[float],
    seed: int,
    threads: int,
) -> List[Tuple[str, float]]:
    obs_by_candidate = [obs[cb] for cb in candidates]
    threads = max(1, min(threads, sim_n))
    chunks = []
    base = sim_n // threads
    rem = sim_n % threads
    start = 0
    for worker in range(threads):
        n = base + (1 if worker < rem else 0)
        if n > 0:
            chunks.append((start, n, seed, probs, log_probs, depths, obs_by_candidate))
            start += n

    if threads == 1:
        partials = [simulate_lower_counts_chunk(chunks[0])]
    else:
        with mp.Pool(processes=threads) as pool:
            partials = pool.map(simulate_lower_counts_chunk, chunks)

    lower = [0] * len(candidates)
    for part in partials:
        for idx, value in enumerate(part):
            lower[idx] += value

    return [(cb, (1.0 + lower[idx]) / (1.0 + sim_n)) for idx, cb in enumerate(candidates)]


def simulate_rrna_twosided_chunk(
    args: Tuple[int, int, int, float, Sequence[int], Sequence[int]]
) -> List[int]:
    start_sim, n_sims, seed, ambient_rrna_p, total_depths, observed_rrna = args
    max_depth = max(total_depths)
    depth_to_candidate_indices: Dict[int, List[int]] = defaultdict(list)
    for idx, depth in enumerate(total_depths):
        depth_to_candidate_indices[depth].append(idx)

    observed_delta = [
        abs((observed_rrna[idx] / total_depths[idx]) - ambient_rrna_p)
        for idx in range(len(total_depths))
    ]
    extreme = [0] * len(total_depths)
    for isim in range(start_sim, start_sim + n_sims):
        rng = random.Random(seed * (isim + 1) + 104729)
        rrna_count = 0
        for depth in range(1, max_depth + 1):
            if rng.random() < ambient_rrna_p:
                rrna_count += 1
            for cand_idx in depth_to_candidate_indices.get(depth, ()):
                sim_delta = abs((rrna_count / depth) - ambient_rrna_p)
                if sim_delta >= observed_delta[cand_idx]:
                    extreme[cand_idx] += 1
    return extreme


def simulate_rrna_twosided_pvalues(
    candidates: List[str],
    total_depths: List[int],
    observed_rrna: List[int],
    ambient_rrna_p: float,
    sim_n: int,
    seed: int,
    threads: int,
) -> List[Tuple[str, float]]:
    threads = max(1, min(threads, sim_n))
    chunks = []
    base = sim_n // threads
    rem = sim_n % threads
    start = 0
    for worker in range(threads):
        n = base + (1 if worker < rem else 0)
        if n > 0:
            chunks.append((start, n, seed, ambient_rrna_p, total_depths, observed_rrna))
            start += n

    if threads == 1:
        partials = [simulate_rrna_twosided_chunk(chunks[0])]
    else:
        with mp.Pool(processes=threads) as pool:
            partials = pool.map(simulate_rrna_twosided_chunk, chunks)

    extreme = [0] * len(candidates)
    for part in partials:
        for idx, value in enumerate(part):
            extreme[idx] += value

    return [(cb, (1.0 + extreme[idx]) / (1.0 + sim_n)) for idx, cb in enumerate(candidates)]


def fisher_combine(p_list: List[float]) -> float:
    """Standard Fisher's method combining all p-values regardless of value."""
    active = []
    for p in p_list:
        if p <= 0:
            p = 1e-300
        active.append(p)
    stat = -2.0 * sum(math.log(p) for p in active)
    df = 2 * len(active)
    try:
        from scipy.stats import chi2
        return float(chi2.sf(stat, df))
    except ImportError:
        raise SystemExit(
            "Fisher combination requires scipy. Install with: pip install scipy"
        ) from None


def soft_umi_pvalue_thre(umis: int, umi_high: int, umi_low: int, thre_high: float, thre_low: float) -> float:
    """Log-linear interpolation of p-value threshold based on UMI depth (soft-end mode).

    For umis >= umi_high: returns thre_high (lenient).
    For umis <= umi_low:  returns thre_low  (strict).
    For umi_low < umis < umi_high: log-linear interpolation.
    """
    if umi_low <= 0:
        return thre_high
    if umis >= umi_high:
        return thre_high
    if umis <= umi_low:
        return thre_low
    # log-linear interpolation
    log_high = math.log(thre_high)
    log_low = math.log(thre_low)
    frac = (umis - umi_low) / (umi_high - umi_low)
    return math.exp(log_low + frac * (log_high - log_low))


def simulate_lower_counts_gene_to_ko_chunk(
    args: Tuple[int, int, int, Sequence[float], Sequence[float], Dict[int, List[int]], Sequence[int], Sequence[float]]
) -> List[int]:
    """Split-mode gene-level simulation with correct multinomial N! accounting.

    Each draw: pick a gene, then increment ALL KOs mapped to that gene.
    logP increment:
      delta_logP = sum_{i=N_old+1}^{N_new} log(i) + sum_ko [log(p_ko) - log(count_ko)]
    where N_new = N_old + m (m = number of KOs for the drawn gene).

    Precomputed log_factorial is passed for efficient N! difference computation.
    """
    start_sim, n_sims, seed, gene_probs, ko_log_probs, gene_to_ko_indices, depths, obs_by_candidate = args
    max_depth = max(depths)
    cumulative = []
    total = 0.0
    for p in gene_probs:
        total += p
        cumulative.append(total)

    depth_to_candidate_indices: Dict[int, List[int]] = defaultdict(list)
    for idx, depth in enumerate(depths):
        depth_to_candidate_indices[depth].append(idx)

    n_kos = len(ko_log_probs)

    max_kos_per_gene = max((len(v) for v in gene_to_ko_indices.values()), default=1)
    max_n = max_depth * max_kos_per_gene
    log_fact = make_log_factorial(max_n)

    lower = [0] * len(depths)
    for isim in range(start_sim, start_sim + n_sims):
        rng = random.Random(seed * (isim + 1))
        curr_counts = [0] * n_kos
        n_total = 0
        logp = 0.0
        for depth in range(1, max_depth + 1):
            g_idx = bisect.bisect_left(cumulative, rng.random() * total)
            ko_indices = gene_to_ko_indices.get(g_idx)
            if not ko_indices:
                continue

            n_old = n_total
            for ko_idx in ko_indices:
                curr_counts[ko_idx] += 1
                logp += ko_log_probs[ko_idx] - math.log(curr_counts[ko_idx])
            n_new = n_total + len(ko_indices)

            logp += log_fact[n_new] - log_fact[n_old]
            n_total = n_new

            for cand_idx in depth_to_candidate_indices.get(depth, ()):
                if logp < obs_by_candidate[cand_idx]:
                    lower[cand_idx] += 1
    return lower


def simulate_pvalues_gene_to_ko(
    candidates: List[str],
    depths: List[int],
    obs: Dict[str, float],
    sim_n: int,
    gene_probs: List[float],
    ko_log_probs: List[float],
    gene_to_ko_indices: Dict[int, List[int]],
    seed: int,
    threads: int,
) -> List[Tuple[str, float]]:
    obs_by_candidate = [obs[cb] for cb in candidates]
    threads = max(1, min(threads, sim_n))
    chunks = []
    base = sim_n // threads
    rem = sim_n % threads
    start = 0
    for worker in range(threads):
        n = base + (1 if worker < rem else 0)
        if n > 0:
            chunks.append((start, n, seed, gene_probs, ko_log_probs, gene_to_ko_indices, depths, obs_by_candidate))
            start += n

    if threads == 1:
        partials = [simulate_lower_counts_gene_to_ko_chunk(chunks[0])]
    else:
        with mp.Pool(processes=threads) as pool:
            partials = pool.map(simulate_lower_counts_gene_to_ko_chunk, chunks)

    lower = [0] * len(candidates)
    for part in partials:
        for idx, value in enumerate(part):
            lower[idx] += value

    return [(cb, (1.0 + lower[idx]) / (1.0 + sim_n)) for idx, cb in enumerate(candidates)]


# ---------------------------------------------------------------------------
# 主测试流程（含 UMI 上下界过滤 + 动态 p 阈值 + no-bh + umi_thre retain）
# ---------------------------------------------------------------------------

def test_kegg_ko(
    cell_gene_counts: CellFeatureCounts,
    umi_per_cb: Counter,
    gene_to_ko: Dict[str, List[str]],
    rrna_genes: Optional[set],
    pvalue_thre: float,
    sim_n: int,
    ambient_pseudocount: float,
    seed: int,
    threads: int,
    ind_min: int,
    ind_max: int,
    umi_lower: int,
    umi_upper: int,
    target_barcodes: Optional[set],
    min_non_rrna_genes: int = 0,
    soft_umi_thre: bool = False,
    soft_umi_high: int = 1000,
    soft_umi_low: int = 100,
    soft_thre_high: float = 0.01,
    soft_thre_low: float = 1e-5,
    umi_thre: Optional[int] = None,
    multi_ko: str = "split",
    no_bh: bool = False,
    omics_mode: bool = False,
    log_likelihood: bool = False,
) -> Tuple[List[Dict[str, object]], Counter, Dict[str, object]]:
    print("[ko] converting non-rRNA gene counts to KO counts", file=sys.stderr, flush=True)
    cell_ko_counts, cell_non_rrna_counts, ko_stats = gene_counts_to_ko_counts(cell_gene_counts, gene_to_ko, rrna_genes)

    print("[test_ko] sorting barcodes", file=sys.stderr, flush=True)
    sorted_cb = sorted(umi_per_cb, key=lambda cb: (-umi_per_cb[cb], cb))
    n_cb = len(sorted_cb)
    rows: Dict[str, Dict[str, object]] = {}
    for rank, cb in enumerate(sorted_cb, 1):
        ko_counts = cell_ko_counts.get(cb, Counter())
        non_rrna_counts = cell_non_rrna_counts.get(cb, Counter())
        rrna_umis = sum(count for gene, count in cell_gene_counts[cb].items() if gene in rrna_genes)
        total_umis = umi_per_cb[cb]
        rows[cb] = {
            "barcode": cb,
            "rank": rank,
            "umis": total_umis,
            "rrna_umis": rrna_umis,
            "rrna_pct": f"{(rrna_umis / total_umis * 100.0) if total_umis else 0.0:.6g}",
            "ko_umis": sum(ko_counts.values()),
            "non_rrna_umis": sum(non_rrna_counts.values()),
            "non_rrna_genes": len(non_rrna_counts),
            "genes": len(cell_gene_counts[cb]),
            "ko_features": len(ko_counts),
            "call": "empty",
            "reason": "not_tested",
            "pvalue_rrnapct": "",
            "padj_rrnapct": "",
            "pvalue_ko": "",
            "padj_ko": "",
            "pvalue_non_rrna": "",
            "padj_non_rrna": "",
            "pvalue_call": "",
            "padj_call": "",
            "pvalue_combined": "",
            "padj_combined": "",
            "pvalue_thre_used": "",
        }

    meta: Dict[str, object] = {
        "n_barcodes": n_cb,
        "mode": "test_kegg_ko",
        "pvalue_thre": pvalue_thre,
        "umi_lower": umi_lower,
        "umi_upper": umi_upper,
    }
    if n_cb <= ind_min:
        meta["warning"] = f"n_barcodes={n_cb} <= ind_min={ind_min}; cannot build ambient profile"
        return [rows[cb] for cb in sorted_cb], ko_stats, meta

    ambient_end = min(n_cb, ind_max)
    ambient_sources = sorted_cb[ind_min:ambient_end]
    print(
        f"[test_ko] building ambient profiles from ranks {ind_min}..{ambient_end} ({len(ambient_sources)} barcodes)",
        file=sys.stderr,
        flush=True,
    )

    # Ambient KO profile
    ambient_counts = Counter()
    all_kos = set()
    for cb, counts in cell_ko_counts.items():
        all_kos.update(counts)
    for cb in ambient_sources:
        ambient_counts.update(cell_ko_counts.get(cb, Counter()))

    if not ambient_counts:
        meta["warning"] = "ambient barcode range has no KO counts; tests skipped"
        return [rows[cb] for cb in sorted_cb], ko_stats, meta

    kos = sorted(all_kos)
    total_ambient_ko = sum(ambient_counts.values()) + ambient_pseudocount * len(kos)
    ko_probs = [(ambient_counts[ko] + ambient_pseudocount) / total_ambient_ko for ko in kos]
    ko_log_probs_list = [math.log(p) for p in ko_probs]
    ko_log_probs = dict(zip(kos, ko_log_probs_list))
    ko_floor_logp = math.log(ambient_pseudocount / total_ambient_ko)

    # Ambient gene profile for third test — controlled by omics_mode
    test_label = "total-RNA" if omics_mode else "non-rRNA"
    ambient_test_counts = Counter()
    for cb in ambient_sources:
        if omics_mode:
            ambient_test_counts.update(cell_gene_counts.get(cb, Counter()))
        else:
            ambient_test_counts.update(cell_non_rrna_counts.get(cb, Counter()))
    all_test_genes = set(ambient_test_counts)

    has_third_test = len(all_test_genes) > 0 and sum(ambient_test_counts.values()) > 0
    if has_third_test:
        test_genes_list = sorted(all_test_genes)
        total_ambient_test = sum(ambient_test_counts.values()) + ambient_pseudocount * len(test_genes_list)
        test_probs = [(ambient_test_counts[g] + ambient_pseudocount) / total_ambient_test for g in test_genes_list]
        test_log_probs_list = [math.log(p) for p in test_probs]
        test_log_probs = dict(zip(test_genes_list, test_log_probs_list))
        test_floor_logp = math.log(ambient_pseudocount / total_ambient_test)
        meta[f"ambient_{test_label}_genes"] = len(test_genes_list)
        meta[f"ambient_{test_label}_umis"] = sum(ambient_test_counts.values())

    # Ambient total-gene rRNA percentage
    ambient_total = 0
    ambient_rrna = 0
    for cb in ambient_sources:
        counts = cell_gene_counts[cb]
        ambient_total += sum(counts.values())
        ambient_rrna += sum(count for gene, count in counts.items() if gene in rrna_genes)
    if ambient_total <= 0:
        meta["warning"] = "ambient barcode range has no total gene counts; tests skipped"
        return [rows[cb] for cb in sorted_cb], ko_stats, meta
    ambient_rrna_p = ambient_rrna / ambient_total
    meta["ambient_rrna_pct"] = ambient_rrna_p * 100.0
    meta["ambient_rrna_umis"] = ambient_rrna
    meta["ambient_total_umis"] = ambient_total

    test_barcodes = sorted_cb
    if target_barcodes is not None:
        test_barcodes = [cb for cb in sorted_cb if cb in target_barcodes]
        print(f"[test_ko] restricting test to {len(test_barcodes)} barcodes from -i BAM", file=sys.stderr, flush=True)

    # -----------------------------------------------------------------------
    # UMI 上下界快速过滤
    # -----------------------------------------------------------------------
    above_upper: List[str] = []
    below_lower: List[str] = []
    middle_barcodes: List[str] = []
    for cb in test_barcodes:
        umis = int(rows[cb]["umis"])
        if umis >= umi_upper:
            above_upper.append(cb)
            # 不提前设为 cell，参与 p 计算，padj 后强制为 0
        elif umis <= umi_lower:
            below_lower.append(cb)
            rows[cb]["call"] = "empty"
            rows[cb]["reason"] = "below_lower_bound"
        else:
            middle_barcodes.append(cb)

    print(
        f"[test_ko] UMI bounds: above_upper({umi_upper}+)={len(above_upper):,}, "
        f"below_lower(<={umi_lower})={len(below_lower):,}, "
        f"middle={len(middle_barcodes):,}",
        file=sys.stderr,
        flush=True,
    )

    # non-rRNA 基因数过滤
    dropped_low_genes: List[str] = []
    if min_non_rrna_genes > 0:
        for cb in list(above_upper):
            if int(rows[cb]["non_rrna_genes"]) < min_non_rrna_genes:
                dropped_low_genes.append(cb)
                rows[cb]["call"] = "empty"
                rows[cb]["reason"] = "below_min_non_rrna_genes"
                above_upper.remove(cb)
        for cb in list(middle_barcodes):
            if int(rows[cb]["non_rrna_genes"]) < min_non_rrna_genes:
                dropped_low_genes.append(cb)
                rows[cb]["call"] = "empty"
                rows[cb]["reason"] = "below_min_non_rrna_genes"
                middle_barcodes.remove(cb)
        print(
            f"[test_ko] min_non_rrna_genes={min_non_rrna_genes}: dropped {len(dropped_low_genes):,} barcodes",
            file=sys.stderr,
            flush=True,
        )

    # above_upper 的 barcode 也参与 p 计算（但不计入 candidates 条件过滤）
    # 将 above_upper 和 middle 合并形成完整的 tested 集合
    tested_barcodes = above_upper + middle_barcodes
    candidates = [cb for cb in tested_barcodes if int(rows[cb]["ko_umis"]) > 0 or int(rows[cb]["non_rrna_umis"]) > 0]
    meta["n_candidates"] = len(candidates)
    meta["n_tested_barcodes"] = len(test_barcodes)
    meta["ambient_rank_start_0based"] = ind_min
    meta["ambient_rank_end_0based"] = ambient_end
    meta["ambient_barcodes"] = len(ambient_sources)
    meta["ambient_ko_features"] = len(kos)
    meta["ambient_ko_umis"] = sum(ambient_counts.values())
    meta["n_above_upper"] = len(above_upper)
    meta["n_below_lower"] = len(below_lower)
    if min_non_rrna_genes > 0:
        meta["min_non_rrna_genes"] = min_non_rrna_genes
        meta["n_dropped_low_genes"] = len(dropped_low_genes)
    if not candidates:
        return [rows[cb] for cb in sorted_cb], ko_stats, meta

    # rRNA two-sided p-value simulation
    print(
        f"[test_ko] simulating two-sided rRNA pct test for {len(candidates)} barcodes "
        f"(ambient rRNA pct={ambient_rrna_p * 100.0:.4g}, sim_n={sim_n}, threads={threads})",
        file=sys.stderr,
        flush=True,
    )
    rrna_pvals = simulate_rrna_twosided_pvalues(
        candidates,
        [umi_per_cb[cb] for cb in candidates],
        [int(rows[cb]["rrna_umis"]) for cb in candidates],
        ambient_rrna_p,
        sim_n,
        seed,
        threads,
    )

    # KO EmptyDrops simulation
    ko_candidates = [cb for cb in candidates if int(rows[cb]["ko_umis"]) > 0]
    if ko_candidates:
        if multi_ko == "split":
            # Gene-level simulation for split mode
            ambient_gene_counts = Counter()
            for cb in ambient_sources:
                for gene, count in cell_gene_counts[cb].items():
                    if gene in rrna_genes:
                        continue
                    if gene in gene_to_ko:
                        ambient_gene_counts[gene] += count

            genes_with_ko = sorted(ambient_gene_counts.keys())
            if not genes_with_ko:
                ko_pvals = []
            else:
                total_ambient_gene = sum(ambient_gene_counts.values()) + ambient_pseudocount * len(genes_with_ko)
                gene_probs = [(ambient_gene_counts[g] + ambient_pseudocount) / total_ambient_gene for g in genes_with_ko]

                ko_to_idx = {ko: i for i, ko in enumerate(kos)}
                gene_to_ko_indices: Dict[int, List[int]] = {}
                for g_idx, gene in enumerate(genes_with_ko):
                    indices = []
                    for ko in gene_to_ko.get(gene, []):
                        if ko in ko_to_idx:
                            indices.append(ko_to_idx[ko])
                    if indices:
                        gene_to_ko_indices[g_idx] = indices

                split_depths = []
                for cb in ko_candidates:
                    mol_count = sum(
                        v for g, v in cell_gene_counts.get(cb, Counter()).items()
                        if g not in (rrna_genes or set()) and g in gene_to_ko
                    )
                    split_depths.append(mol_count)

                print(
                    f"[test_ko] split mode: gene-level simulation with {len(genes_with_ko)} genes, "
                    f"max mol depth={max(split_depths) if split_depths else 0}",
                    file=sys.stderr,
                    flush=True,
                )

                # Compute obs_ko for split mode (same multinomial logP as non-split)
                max_ko = max(int(rows[cb]["ko_umis"]) for cb in ko_candidates)
                log_factorial_ko = make_log_factorial(max_ko)
                obs_ko = {
                    cb: log_multinomial(cell_ko_counts[cb], ko_log_probs, log_factorial_ko, ko_floor_logp)
                    for cb in ko_candidates
                }
                if log_likelihood:
                    for cb, val in obs_ko.items():
                        rows[cb]["_loglik_ko"] = val

                ko_pvals = simulate_pvalues_gene_to_ko(
                    ko_candidates,
                    split_depths,
                    obs_ko,
                    sim_n,
                    gene_probs,
                    ko_log_probs_list,
                    gene_to_ko_indices,
                    seed,
                    threads,
                )
        else:
            max_ko = max(int(rows[cb]["ko_umis"]) for cb in ko_candidates)
            print(
                f"[test_ko] simulating KO EmptyDrops for {len(ko_candidates)} barcodes "
                f"(max KO UMI={max_ko}, sim_n={sim_n}, threads={threads})",
                file=sys.stderr,
                flush=True,
            )
            log_factorial_ko = make_log_factorial(max_ko)
            obs_ko = {
                cb: log_multinomial(cell_ko_counts[cb], ko_log_probs, log_factorial_ko, ko_floor_logp)
                for cb in ko_candidates
            }
            if log_likelihood:
                for cb, val in obs_ko.items():
                    rows[cb]["_loglik_ko"] = val
            ko_pvals = simulate_pvalues(
                ko_candidates,
                [int(rows[cb]["ko_umis"]) for cb in ko_candidates],
                obs_ko,
                sim_n,
                ko_probs,
                ko_log_probs_list,
                seed,
                threads,
            )
    else:
        ko_pvals = []

    # Gene-level EmptyDrops simulation (third test) — controlled by omics_mode
    if omics_mode:
        third_candidates = [cb for cb in candidates if int(rows[cb]["umis"]) > 0]
        third_test_counts = cell_gene_counts
        third_depth_key = "umis"
    else:
        third_candidates = [cb for cb in candidates if int(rows[cb]["non_rrna_umis"]) > 0]
        third_test_counts = cell_non_rrna_counts
        third_depth_key = "non_rrna_umis"
    if third_candidates and has_third_test:
        max_um = max(int(rows[cb][third_depth_key]) for cb in third_candidates)
        print(
            f"[test_ko] simulating {test_label} gene EmptyDrops for {len(third_candidates)} barcodes "
            f"(max {test_label} UMI={max_um}, sim_n={sim_n}, threads={threads})",
            file=sys.stderr,
            flush=True,
        )
        log_factorial_um = make_log_factorial(max_um)
        obs_um = {
            cb: log_multinomial(third_test_counts[cb], test_log_probs, log_factorial_um, test_floor_logp)
            for cb in third_candidates
        }
        if log_likelihood:
            for cb, val in obs_um.items():
                rows[cb]["_loglik_gene"] = val
        um_pvals = simulate_pvalues(
            third_candidates,
            [int(rows[cb][third_depth_key]) for cb in third_candidates],
            obs_um,
            sim_n,
            test_probs,
            test_log_probs_list,
            seed,
            threads,
        )
    else:
        um_pvals = []

    # Build p-value dicts (default p=1.0 for barcodes not in a sub-test)
    rrna_p = dict(rrna_pvals)
    ko_p = dict(ko_pvals)
    um_p = dict(um_pvals)

    rrna_padj = bh_adjust(rrna_pvals)
    ko_padj = bh_adjust(ko_pvals)
    um_padj = bh_adjust(um_pvals)

    # Fisher combination: always combine all three p-values (no supportive testing)
    combined_pairs = []
    for cb in candidates:
        p_ko_val = ko_p.get(cb, 1.0)
        p_rrna_val = rrna_p.get(cb, 1.0)
        p_test_val = um_p.get(cb, 1.0)
        combined_p = fisher_combine([p_ko_val, p_rrna_val, p_test_val])
        combined_pairs.append((cb, combined_p))

    # UMI retain / above_upper 逻辑：这些 barcode padj 强制为 0，参与 BH 校正
    retained_bcs = set(above_upper)
    if umi_thre is not None:
        retained_bcs |= {cb for cb in candidates if umi_per_cb[cb] > umi_thre}
    if retained_bcs:
        above_from_upper = len(above_upper)
        above_from_thre = len(retained_bcs) - above_from_upper
        parts = []
        if above_from_upper:
            parts.append(f"UMI>={umi_upper}: {above_from_upper:,}")
        if above_from_thre:
            parts.append(f"UMI>{umi_thre}: {above_from_thre:,}")
        print(
            f"[test_ko] retaining {len(retained_bcs):,} barcodes ({', '.join(parts)}) p=0, participates in BH",
            file=sys.stderr,
            flush=True,
        )

    combined_pairs_all = []
    for cb, p in combined_pairs:
        if cb in retained_bcs:
            combined_pairs_all.append((cb, 0.0))
        else:
            combined_pairs_all.append((cb, p))
    combined_adj = bh_adjust(combined_pairs_all)

    n_cells = 0
    for cb in candidates:
        p_ko_val = ko_p.get(cb, 1.0)
        p_rrna_val = rrna_p.get(cb, 1.0)
        p_test_val = um_p.get(cb, 1.0)
        call_pvalue = fisher_combine([p_ko_val, p_rrna_val, p_test_val])
        call_padj = combined_adj[cb]
        umis = int(rows[cb]["umis"])
        if soft_umi_thre:
            effective_thre = soft_umi_pvalue_thre(umis, soft_umi_high, soft_umi_low, soft_thre_high, soft_thre_low)
        else:
            effective_thre = pvalue_thre
        rows[cb]["pvalue_rrnapct"] = f"{p_rrna_val:.8g}"
        rows[cb]["padj_rrnapct"] = f"{rrna_padj.get(cb, 1.0):.8g}"
        rows[cb]["pvalue_ko"] = f"{p_ko_val:.8g}"
        rows[cb]["padj_ko"] = f"{ko_padj.get(cb, 1.0):.8g}"
        rows[cb]["pvalue_non_rrna"] = f"{p_test_val:.8g}"
        rows[cb]["padj_non_rrna"] = f"{um_padj.get(cb, 1.0):.8g}"
        rows[cb]["pvalue_call"] = f"{call_pvalue:.8g}"
        rows[cb]["padj_call"] = f"{call_padj:.8g}"
        rows[cb]["pvalue_combined"] = f"{call_pvalue:.8g}"
        rows[cb]["padj_combined"] = f"{call_padj:.8g}"
        rows[cb]["pvalue_thre_used"] = f"{effective_thre:.8g}"
        rows[cb]["reason"] = "tested"

        if cb in above_upper and cb in retained_bcs:
            rows[cb]["call"] = "cell"
            rows[cb]["reason"] = "above_upper_bound"
            n_cells += 1
        elif cb in retained_bcs:
            rows[cb]["call"] = "cell"
            rows[cb]["reason"] = "umi_retained"
            n_cells += 1
        elif no_bh:
            if call_pvalue < effective_thre:
                rows[cb]["call"] = "cell"
                rows[cb]["reason"] = "combined_pvalue_no_bh"
                n_cells += 1
        elif call_padj < effective_thre:
            rows[cb]["call"] = "cell"
            rows[cb]["reason"] = "combined_padj"
            n_cells += 1

    meta["n_cells_total"] = sum(1 for cb in sorted_cb if rows[cb]["call"] == "cell")
    meta["simulation_threads"] = threads
    meta["omics_mode"] = omics_mode
    meta["test_label"] = test_label
    bh_label = "no_bh" if no_bh else "bh"
    if soft_umi_thre:
        meta["call_rule"] = (
            f"soft_umi_thre(umi_high={soft_umi_high},umi_low={soft_umi_low},"
            f"thre_high={soft_thre_high},thre_low={soft_thre_low}); "
            f"fisher_combined_{bh_label}_pvalue<{pvalue_thre}; "
            f"fisher_combine(p_ko+p_rrna+p_{test_label} unconditional); "
            f"umi>={umi_upper} auto-cell, umi<={umi_lower} auto-empty"
        )
    else:
        meta["call_rule"] = (
            f"fisher_combined_{bh_label}_pvalue<{pvalue_thre}; "
            f"fisher_combine(p_ko+p_rrna+p_{test_label} unconditional); "
            f"umi>={umi_upper} auto-cell, umi<={umi_lower} auto-empty"
        )
    if umi_thre is not None:
        meta["call_rule"] += f"; umi>{umi_thre} forced retained"
    return [rows[cb] for cb in sorted_cb], ko_stats, meta


# ---------------------------------------------------------------------------
# 输出
# ---------------------------------------------------------------------------

def write_tsv(path: str, rows, fields: List[str]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_outputs(out_prefix: str, rows: List[Dict[str, object]], stats: Counter, meta: Dict[str, object], log_likelihood: bool = False) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(out_prefix)) or ".", exist_ok=True)
    fields = [
        "barcode",
        "rank",
        "umis",
        "rrna_umis",
        "rrna_pct",
        "ko_umis",
        "non_rrna_umis",
        "non_rrna_genes",
        "genes",
        "ko_features",
        "call",
        "reason",
        "pvalue_rrnapct",
        "padj_rrnapct",
        "pvalue_ko",
        "padj_ko",
        "pvalue_non_rrna",
        "padj_non_rrna",
        "pvalue_call",
        "padj_call",
        "pvalue_combined",
        "padj_combined",
        "pvalue_thre_used",
    ]
    tested_rows = [row for row in rows if row.get("reason") in ("tested", "combined_padj", "combined_pvalue_no_bh", "umi_retained", "above_upper_bound", "below_min_non_rrna_genes")]
    write_tsv(out_prefix + ".barcodes.tsv", tested_rows, fields)

    cells = [str(row["barcode"]) for row in rows if row["call"] == "cell"]
    with open(out_prefix + ".cells.txt", "w", encoding="utf-8") as handle:
        for cb in cells:
            handle.write(cb + "\n")

    with open(out_prefix + ".summary.txt", "w", encoding="utf-8") as handle:
        for key in sorted(stats):
            handle.write(f"{key}\t{stats[key]}\n")
        for key in sorted(meta):
            handle.write(f"{key}\t{meta[key]}\n")
        handle.write(f"cells_called\t{len(cells)}\n")

    if log_likelihood:
        test_label = meta.get("test_label", "non-rRNA")
        loglik_fields = [
            "barcode",
            "rank",
            "umis",
            "ko_umis",
            "loglik_ko",
            f"loglik_{test_label}",
        ]
        loglik_rows = []
        for row in tested_rows:
            loglik_rows.append({
                "barcode": row["barcode"],
                "rank": row["rank"],
                "umis": row["umis"],
                "ko_umis": row["ko_umis"],
                "loglik_ko": f"{row.get('_loglik_ko', ''):.8g}" if row.get("_loglik_ko") is not None else "",
                f"loglik_{test_label}": f"{row.get('_loglik_gene', ''):.8g}" if row.get("_loglik_gene") is not None else "",
            })
        write_tsv(out_prefix + ".loglik.tsv", loglik_rows, loglik_fields)
        print(f"[output] wrote log-likelihood to {out_prefix}.loglik.tsv", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# 主函数
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="MVAcaller (multivariate ambient caller): KEGG-KO / EmptyDrops-style cell calling for microbial scRNA-seq. In-memory and SQLite modes, with UMI upper/lower bound filters and a default soft UMI-aware P-threshold."
    )
    # 输入
    parser.add_argument("-i", "--input", default=None,
                        help="Input gene-tagged BAM/SAM (optional when an existing --build-sqlite database is reused)")
    parser.add_argument("-g", "--gff", required=True, help="GFF3 annotation (used to identify rRNA genes)")
    parser.add_argument("--kegg-map", required=True, help="Gene -> KEGG KO mapping TSV (two columns: gene_id\tKO)")

    # 输出
    parser.add_argument("-o", "--out-prefix", required=True, help="Output prefix (produces <prefix>.barcodes.tsv, <prefix>.cells.txt, <prefix>.summary.txt)")
    parser.add_argument(
        "--build-sqlite",
        default=None,
        metavar="PATH",
        help="SQLite mode: build / load molecule DB at this path; omit for in-memory mode",
    )
    parser.add_argument(
        "--reuse-sqlite",
        action="store_true",
        help="Reuse an existing --build-sqlite file instead of rebuilding",
    )

    # BAM 解析参数
    parser.add_argument(
        "--qname-regex",
        default=r"_(?P<cb>[ACGTN]+)_(?P<umi>[ACGTN]+)$",
        help="Regular expression for extracting CB and UMI from qname, with named groups (?P<cb>...) (?P<umi>...)",
    )
    parser.add_argument("--use-cb-tag", action="store_true", help="Extract CB/UMI from SAM tags (CB:Z: + UB:Z: or UR:Z:) instead of qname")
    parser.add_argument("--qname-simple", action=argparse.BooleanOptionalAction, default=True, help="Use the simple qname regex _([^_]+)_(.+)$ (default ON; turn off with --no-qname-simple to use --qname-regex)")
    parser.add_argument("--gene-delimiter", default="#", help="Delimiter between contig and gene id in RNAME; default: #")
    parser.add_argument("--min-mapq", type=int, default=0, help="Minimum MAPQ; default: 0")
    parser.add_argument("--keep-secondary", action="store_true", help="Keep secondary / supplementary alignments")
    parser.add_argument("--batch-size", type=int, default=200000, help="SQLite insert batch size; default: 200000")
    parser.add_argument("--progress-every", type=int, default=1000000, help="Print progress every N reads; 0 disables")

    # EmptyDrops 参数
    parser.add_argument("--ind-min", type=int, default=None, help="Start rank (inclusive) of the ambient-profile barcode range. Default: auto = Q3 rank (bottom 25%% of barcodes by UMI)")
    parser.add_argument("--ind-max", type=int, default=570000000, help="End rank (inclusive) of the ambient-profile barcode range; default 570000000 (effectively up to the end of the tail)")
    parser.add_argument("--sim-n", type=int, default=10000, help="Number of multinomial simulations per barcode")
    parser.add_argument("--threads", type=int, default=1, help="Worker threads for simulations")
    parser.add_argument("--ambient-pseudocount", type=float, default=1.0, help="Ambient-profile pseudocount; default: 1.0")
    parser.add_argument("--seed", type=int, default=19760110, help="Random seed")
    parser.add_argument("--pvalue-thre", type=float, default=0.01, help="Combined P-value threshold (post-BH) for cell calls; default: 0.01")
    parser.add_argument("--no-bh", action="store_true", help="Disable Benjamini-Hochberg correction and call cells directly on the raw Fisher combined P-value")
    parser.add_argument(
        "--umi-thre",
        type=int,
        default=None,
        help="UMI retention threshold: barcodes with UMI > this are forced to cell (P=0, included in BH correction)",
    )

    # UMI 上下界过滤（核心新增）
    parser.add_argument(
        "--umi-lower",
        type=int,
        default=None,
        help="UMI lower bound (inclusive). Barcodes with UMI <= this are skipped as empty. Default: 300 (single-species) / 100 (--mode-omics)",
    )
    parser.add_argument(
        "--umi-upper",
        type=int,
        default=10000,
        help="UMI upper bound (inclusive). Barcodes with UMI >= this are forced to cell without computing P. Default 10000",
    )
    parser.add_argument(
        "--min-non-rrna-genes",
        type=int,
        default=0,
        help="Lower bound on the number of detected non-rRNA genes; barcodes below this are skipped and treated as empty. Default 0 (off)",
    )

    # KO 参数
    parser.add_argument(
        "--multi-ko",
        choices=["split", "first", "drop"],
        default="split",
        help="How to handle genes annotated to multiple KOs; default: split (counts split equally)",
    )
    parser.add_argument(
        "--target-barcodes",
        default=None,
        help="Optional file, one barcode per line; restrict testing to these barcodes only",
    )
    parser.add_argument(
        "--soft-umi-thre",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Soft-end P-threshold mode (ON by default); threshold slides along log(UMI), stricter for low-UMI barcodes. Turn off with --no-soft-umi-thre",
    )
    parser.add_argument("--soft-umi-high", type=int, default=1000, help="soft-end high UMI anchor (>= this uses relaxed threshold); default 1000")
    parser.add_argument("--soft-umi-low", type=int, default=100, help="soft-end low UMI anchor (<= this uses strict threshold); default 100")
    parser.add_argument("--soft-thre-high", type=float, default=0.01, help="P-threshold at the high-UMI end; default 0.01")
    parser.add_argument("--soft-thre-low", type=float, default=1e-5, help="P-threshold at the low-UMI end; default 1e-5")
    parser.add_argument(
        "--mode-omics",
        action="store_true",
        default=False,
        help="Omics mode: replace the third (non-rRNA composition) test with a total-RNA (including rRNA) test",
    )
    parser.add_argument(
        "--log-likelihood",
        action="store_true",
        default=False,
        help="Write per-barcode log-multinomial likelihood to <prefix>.loglik.tsv (off by default)",
    )

    args = parser.parse_args()

    # -i is not needed when counts are already available as a reused SQLite database
    if args.reuse_sqlite and not args.build_sqlite:
        parser.error("--reuse-sqlite requires --build-sqlite PATH")
    reusing_sqlite = bool(
        args.build_sqlite and args.reuse_sqlite and os.path.exists(args.build_sqlite)
    )
    if not reusing_sqlite and not args.input:
        parser.error("-i/--input is required unless an existing --build-sqlite database is reused")

    # Mode-dependent default for --umi-lower
    if args.umi_lower is None:
        args.umi_lower = 100 if args.mode_omics else 300
        print(
            f"[defaults] --umi-lower = {args.umi_lower} "
            f"(mode={'omics' if args.mode_omics else 'single-species'})",
            file=sys.stderr, flush=True,
        )

    # UMI 上下界验证
    if args.umi_upper > 0 and args.umi_lower >= args.umi_upper:
        print(f"[error] --umi-lower ({args.umi_lower}) must be < --umi-upper ({args.umi_upper})", file=sys.stderr)
        return 1

    # soft-end 参数验证
    if args.soft_umi_thre:
        if args.soft_umi_low <= 0:
            print("[error] --soft-umi-low must be > 0", file=sys.stderr)
            return 1
        if args.soft_umi_high <= args.soft_umi_low:
            print(f"[error] --soft-umi-high ({args.soft_umi_high}) must be > --soft-umi-low ({args.soft_umi_low})", file=sys.stderr)
            return 1
        if not (0.0 < args.soft_thre_high <= 1.0) or not (0.0 < args.soft_thre_low <= args.soft_thre_high):
            print(f"[error] invalid soft threshold range: --soft-thre-low ({args.soft_thre_low}) must be in (0, --soft-thre-high ({args.soft_thre_high})]", file=sys.stderr)
            return 1
        print(
            f"[soft_umi_thre] enabled: UMI {args.soft_umi_low}->{args.soft_umi_high}, "
            f"thre {args.soft_thre_low}->{args.soft_thre_high}",
            file=sys.stderr, flush=True,
        )

    # 正则 (precedence: --use-cb-tag > --no-qname-simple + --qname-regex > --qname-simple default)
    if args.use_cb_tag:
        regex = None
    elif not args.qname_simple:
        regex = re.compile(args.qname_regex)
    else:
        regex = re.compile(r"_(?P<cb>[^_]+)_(?P<umi>.+)$")

    # -----------------------------------------------------------------------
    # 模式选择：内存 or SQLite
    # -----------------------------------------------------------------------
    if args.build_sqlite:
        if args.reuse_sqlite and os.path.exists(args.build_sqlite):
            print(f"[sqlite] reusing existing {args.build_sqlite}", file=sys.stderr, flush=True)
            cell_gene_counts, umi_per_cb, stats = load_gene_counts_sqlite(args.build_sqlite)
        else:
            cell_gene_counts, umi_per_cb, stats = build_sqlite(
                args.input,
                regex,
                args.gene_delimiter,
                args.min_mapq,
                not args.keep_secondary,
                args.build_sqlite,
                args.batch_size,
                args.progress_every,
                use_cb_tag=args.use_cb_tag,
            )
    else:
        cell_gene_counts, umi_per_cb, stats = build_counts_from_bam(
            args.input,
            regex,
            args.gene_delimiter,
            args.min_mapq,
            not args.keep_secondary,
            args.progress_every,
            use_cb_tag=args.use_cb_tag,
        )

    # 注释
    rrna_genes = load_rrna_genes(args.gff)
    stats["rrna_genes"] = len(rrna_genes)

    gene_to_ko = load_gene_to_ko(args.kegg_map, args.multi_ko)
    stats["kegg_mapped_genes"] = len(gene_to_ko)
    stats["multi_ko_mode"] = args.multi_ko

    # target barcodes
    target_barcodes = None
    if args.target_barcodes:
        target_barcodes = set()
        with open(args.target_barcodes, "r", encoding="utf-8") as fh:
            for line in fh:
                cb = line.strip()
                if cb:
                    target_barcodes.add(cb)
        stats["target_barcodes"] = len(target_barcodes)
        missing = target_barcodes - set(umi_per_cb)
        if missing:
            print(f"[target] {len(missing):,} target barcodes not in input, adding empty entries", file=sys.stderr, flush=True)
            for cb in missing:
                cell_gene_counts[cb] = Counter()
                umi_per_cb[cb] = 0
            stats["target_barcodes_missing_from_db"] = len(missing)

    # Auto-resolve --ind-min: Q3 rank = bottom 25% of barcodes by UMI (the ambient tail)
    if args.ind_min is None:
        n_cb = len(umi_per_cb)
        args.ind_min = max(int(n_cb * 0.75), 1)
        print(
            f"[defaults] --ind-min = {args.ind_min} (auto Q3 rank, "
            f"n_barcodes={n_cb}, bottom 25% tail)",
            file=sys.stderr, flush=True,
        )

    # 测试
    rows, ko_stats, meta = test_kegg_ko(
        cell_gene_counts,
        umi_per_cb,
        gene_to_ko,
        rrna_genes,
        args.pvalue_thre,
        args.sim_n,
        args.ambient_pseudocount,
        args.seed,
        args.threads,
        args.ind_min,
        args.ind_max,
        args.umi_lower,
        args.umi_upper,
        target_barcodes,
        min_non_rrna_genes=args.min_non_rrna_genes,
        soft_umi_thre=args.soft_umi_thre,
        soft_umi_high=args.soft_umi_high,
        soft_umi_low=args.soft_umi_low,
        soft_thre_high=args.soft_thre_high,
        soft_thre_low=args.soft_thre_low,
        umi_thre=args.umi_thre,
        multi_ko=args.multi_ko,
        no_bh=args.no_bh,
        omics_mode=args.mode_omics,
        log_likelihood=args.log_likelihood,
    )
    stats.update(ko_stats)
    write_outputs(args.out_prefix, rows, stats, meta, log_likelihood=args.log_likelihood)

    print(f"barcodes: {stats['barcodes_detected']}")
    print(f"KO features: {stats['ko_features_detected']}")
    print(f"cells called: {sum(1 for row in rows if row['call'] == 'cell')}")
    if "warning" in meta:
        print(f"warning: {meta['warning']}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
