#!/usr/bin/env python3
"""
Convert UHGG GTF (key "value" attributes) to GFF3-style (key=value),
so run_kegg_ko_emptydropsx.py can parse it.

Key changes:
1. gene_id "xxx" -> ID=xxx  (script requires attrs.get("ID"))
2. key "value"  -> key=value (script uses split("=") parsing)
"""
import re
import sys

GTF_PATH = sys.argv[1] if len(sys.argv) > 1 else "/Users/apple/Downloads/all_species.gtf"
OUT_PATH = sys.argv[2] if len(sys.argv) > 2 else "/Users/apple/Downloads/all_species.gff3"

ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')

with open(GTF_PATH) as fin, open(OUT_PATH, "w") as fout:
    for line in fin:
        if line.startswith("#"):
            fout.write(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            fout.write(line)
            continue

        attrs = dict(ATTR_RE.findall(fields[8]))
        gene_id = attrs.pop("gene_id", "")

        parts = []
        if gene_id:
            parts.append(f"ID={gene_id}")
        for k, v in attrs.items():
            parts.append(f"{k}={v}")
        fields[8] = ";".join(parts)
        fout.write("\t".join(fields) + "\n")
