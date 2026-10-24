# Summarise an MVAcaller *.barcodes.tsv by barcode class.
# Classes follow the benchmark labelling: "EMPTY" prefix = simulated empty droplet,
# "#small" suffix = downsampled cell, everything else = large cell.
# Prints one TSV line:
#   called_combined_large called_combined_small called_combined_empty \
#   called_nonrrna_large  called_nonrrna_small  called_nonrrna_empty \
#   total_large           total_small           total_empty
BEGIN { FS = OFS = "\t" }
NR == 1 {
  for (i = 1; i <= NF; i++) h[$i] = i
  if (!("barcode" in h) || !("call" in h)) {
    print "missing barcode/call column" > "/dev/stderr"; exit 2
  }
  next
}
{
  cb = $h["barcode"]
  cls = (cb ~ /EMPTY/) ? "empty" : ((cb ~ /#small/) ? "small" : "large")
  tot[cls]++
  if ($h["call"] == "cell") comb[cls]++
  # an empty p-value means the barcode was never tested, not p = 0
  if ($h["padj_non_rrna"] != "" && ($h["padj_non_rrna"] + 0) < TH) nr[cls]++
  if ($h["padj_combined"] != "" && ($h["padj_combined"] + 0) < TH) all[cls]++
}
END {
  printf "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
         all["large"], all["small"], all["empty"],
         nr["large"],  nr["small"],  nr["empty"],
         tot["large"], tot["small"], tot["empty"]
}
