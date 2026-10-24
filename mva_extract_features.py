#!/usr/bin/env python3
"""
Extract per-gene sequences from a genome FASTA + annotation (GTF or GFF3).

Produces, in the output directory:
  cds.fa       protein-coding genes (CDS)
  trna.fa      tRNA genes
  rrna.fa      rRNA genes
  ncrna.fa     all remaining non-coding genes (ncRNA, tmRNA, CRISPR, SRP_RNA, ...)
  features.fa  cds.fa + trna.fa + rrna.fa + ncrna.fa, used to build the gene-tagged
               alignment index (bowtie2/STAR) consumed by MVAcaller.py

Headers follow the MVAcaller gene-tagged convention:
    ><seqid>#<gene_id>
so that the aligned RNAME is <contig>#<gene_id> and --gene-delimiter '#' (default)
recovers the bare gene_id.

Usage:
    python mva_extract_features.py -g annotation.gtf -f genome.fna -o ref/
"""

import argparse
import os
import re
import sys

# Non-coding biotypes kept in the reference (mobile elements are skipped)
TRNA_BIOTYPES = {"tRNA", "tRNA_gene"}
RRNA_BIOTYPES = {"rRNA", "rRNA_gene"}
NCRNA_BIOTYPES = {
    "ncRNA", "tmRNA", "CRISPR", "CRISPRdr", "CRISPRspacer",
    "RNase_MRP", "RNase_P", "SRP_RNA", "ncRNA_gene",
}
CDS_BIOTYPES = {"protein_coding", "CDS"}

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def parse_attributes(attr_str):
    """Parse the 9th column. Auto-detects GFF3 (key=value) vs GTF (key "value")."""
    attrs = {}
    if "=" in attr_str and '"' not in attr_str:
        for pair in attr_str.split(";"):
            pair = pair.strip()
            if "=" in pair:
                key, value = pair.split("=", 1)
                attrs[key.strip()] = value.strip()
    else:
        for match in re.finditer(r'(\w+)\s+"([^"]*)"', attr_str):
            attrs[match.group(1)] = match.group(2)
    return attrs


def load_fasta(path):
    seqs = {}
    current_id = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(chunks)
                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current_id is not None:
        seqs[current_id] = "".join(chunks)
    return seqs


def main():
    parser = argparse.ArgumentParser(
        description="Extract CDS / tRNA / rRNA / ncRNA sequences and a combined "
                    "gene-tagged feature FASTA from a genome FASTA + GTF/GFF3."
    )
    parser.add_argument("-g", "--gtf", required=True, help="Input annotation (GTF or GFF3)")
    parser.add_argument("-f", "--fna", required=True, help="Input genome FASTA (.fna)")
    parser.add_argument("-o", "--output", default="ref", help="Output directory (default: ref)")
    parser.add_argument(
        "--keep-mobile", action="store_true",
        help="Also keep mobile elements (insertion_sequence, prophage, plasmid, integron); off by default",
    )
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    print("Loading FASTA...", file=sys.stderr, flush=True)
    seqs = load_fasta(args.fna)
    print(f"  Loaded {len(seqs)} sequences", file=sys.stderr, flush=True)

    paths = {name: os.path.join(args.output, name) for name in
             ("cds.fa", "trna.fa", "rrna.fa", "ncrna.fa", "features.fa")}
    handles = {name: open(path, "w") for name, path in paths.items()}
    counts = {"cds": 0, "trna": 0, "rrna": 0, "ncrna": 0, "mobile": 0, "skipped": 0}

    with open(args.gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, _source, ftype, start, end, _score, strand, _frame, attrs_str = fields

            attrs = parse_attributes(attrs_str)
            biotype = attrs.get("gene_biotype", "") or ftype
            gene_id = attrs.get("gene_id", "") or attrs.get("ID", "") or attrs.get("locus_tag", "")
            if not gene_id:
                counts["skipped"] += 1
                continue

            if seqid not in seqs:
                counts["skipped"] += 1
                continue

            subseq = seqs[seqid][int(start) - 1:int(end)]
            if strand == "-":
                subseq = revcomp(subseq)

            record = f">{seqid}#{gene_id}\n{subseq}\n"

            if biotype in CDS_BIOTYPES:
                bucket, key = "cds.fa", "cds"
            elif biotype in TRNA_BIOTYPES:
                bucket, key = "trna.fa", "trna"
            elif biotype in RRNA_BIOTYPES:
                bucket, key = "rrna.fa", "rrna"
            elif biotype in NCRNA_BIOTYPES:
                bucket, key = "ncrna.fa", "ncrna"
            elif args.keep_mobile:
                bucket, key = "ncrna.fa", "mobile"
            else:
                counts["mobile"] += 1
                continue

            handles[bucket].write(record)
            handles["features.fa"].write(record)
            counts[key] += 1

    for handle in handles.values():
        handle.close()

    print("\nDone!", file=sys.stderr, flush=True)
    for key, label in (("cds", "CDS"), ("trna", "tRNA"), ("rrna", "rRNA"),
                       ("ncrna", "ncRNA")):
        print(f"  {label:<6} {counts[key]:>9,} -> {paths[key + '.fa' if key != 'cds' else 'cds.fa']}",
              file=sys.stderr, flush=True)
    total = counts["cds"] + counts["trna"] + counts["rrna"] + counts["ncrna"]
    print(f"  Total  {total:>9,} -> {paths['features.fa']}", file=sys.stderr, flush=True)
    if counts["skipped"]:
        print(f"  Skipped (seqid not found / no gene id): {counts['skipped']}", file=sys.stderr)
    if counts["mobile"]:
        print(f"  Skipped mobile elements: {counts['mobile']}", file=sys.stderr)


if __name__ == "__main__":
    main()
