#!/usr/bin/env bash
# Demo on public reference sequence only (GRCm39, mouse C57BL/6J).
#
# Needs:
#   samtools, blastn, and a BLAST database built from GRCm39
#   python with biopython (see ../requirements.txt)
#
# Build the database once, if you have not:
#   makeblastdb -in GRCm39.fa -dbtype nucl -out ~/genomes/blastdb/GRCm39
#
# Usage:  bash examples/demo_grcm39.sh [GENOME_FASTA] [BLAST_DB] [OUTDIR]
set -euo pipefail

GENOME="${1:-$HOME/genomes/GRCm39.bgz}"
DB="${2:-$HOME/genomes/blastdb/GRCm39}"
OUT="${3:-./demo-out}"
REGION="14:8440000-8470000"          # 30 kb inside the mouse Atxn7 locus
PY="${PYTHON:-python3}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

for t in samtools blastn; do
  command -v "$t" >/dev/null || { echo "missing: $t" >&2; exit 1; }
done
[ -r "$GENOME" ] || { echo "no genome FASTA at $GENOME" >&2; exit 1; }
ls "$DB".n* >/dev/null 2>&1 || { echo "no BLAST database at $DB" >&2; exit 1; }

mkdir -p "$OUT"
echo "1/4  pulling $REGION from $(basename "$GENOME")"
samtools faidx "$GENOME" "$REGION" > "$OUT/locus.fa"

echo "2/4  repeat-density map (this is the slow step, ~1 min)"
"$PY" "$HERE/lrpcr/repeatmap.py" --fasta "$OUT/locus.fa" --db "$DB" \
      --tile 1000 -o "$OUT/repeat_map.tsv"

echo "3/4  taking three 30-mers straight off the reference as demo primers"
"$PY" - "$OUT" <<'PY'
import sys
out = sys.argv[1]
seq = "".join(l.strip() for l in open(f"{out}/locus.fa") if not l.startswith(">")).upper()
# Named by offset only: nothing is implied about their quality.
with open(f"{out}/primers.fa", "w") as fh:
    for pos in (400, 2400, 3400):
        fh.write(f">demo_at_{pos}\n{seq[pos-1:pos+29]}\n")
PY

echo "4/4  3'-anchor screen"
"$PY" "$HERE/lrpcr/anchors.py" --primers "$OUT/primers.fa" --db "$DB" \
      --chunk 3 -o "$OUT/anchors.tsv"

echo
echo "repeat-density map, tile verdicts:"
awk -F'\t' 'NR>1{c[$4]++} END{for (k in c) printf "  %-8s %d tiles\n", k, c[k]}' \
    "$OUT/repeat_map.tsv"
echo
echo "3'-anchor screen:"
column -t -s$'\t' "$OUT/anchors.tsv" | sed 's/^/  /'
echo
cat <<'TXT'
Read it like this: the two signals are independent. A primer can sit in a tile
the map calls clean and still carry thousands of 3' anchors, and vice versa.
Neither number alone predicts a clean gel; a candidate has to pass both, and
self_products is the one that tracks laddering.

The digest tool needs amplicons rather than a locus, so it is not part of this
demo. See docs/METHOD.md.
TXT
