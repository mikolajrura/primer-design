#!/usr/bin/env python3
"""Repeat-density map of a locus, in tiles, before you design anything.

Designing primers first and discarding the repetitive ones afterwards wastes
most of the work: whole kilobase-scale windows of a mammalian locus sit inside
SINE/LINE families, and every candidate from such a window fails. Tiling the
region first and counting how often each tile recurs in the genome tells you
where it is worth running primer3 at all.

Each tile is BLASTed (megablast) against the genome and its hits of at least
--min-hit-len are counted. A unique tile scores 1 (itself). A tile inside a B1
element scores tens of thousands.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def read_single_fasta(path: Path) -> tuple[str, str]:
    name, chunks = None, []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                break  # only the first record
            name = line[1:].split()[0]
        elif name is not None:
            chunks.append(line)
    if name is None:
        raise SystemExit(f"{path}: no FASTA record")
    return name, "".join(chunks).upper()


def parse_region(text: str, seq_len: int) -> tuple[int, int]:
    if not text:
        return 1, seq_len
    try:
        lo_s, hi_s = text.split("-")
        lo, hi = int(lo_s.replace(",", "")), int(hi_s.replace(",", ""))
    except ValueError:
        raise SystemExit(f"--region must look like 99000-130000, got {text!r}")
    if not (1 <= lo < hi <= seq_len):
        raise SystemExit(f"--region {lo}-{hi} outside 1-{seq_len}")
    return lo, hi


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Tile a locus and count how often each tile recurs in a genome.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--fasta", type=Path, required=True,
                    help="locus sequence; only the first record is used")
    ap.add_argument("--db", required=True, help="BLAST nucleotide database (makeblastdb prefix)")
    ap.add_argument("--region", default="", help="1-based inclusive span, e.g. 99000-130000")
    ap.add_argument("--tile", type=int, default=500, help="tile size in bp")
    ap.add_argument("--min-hit-len", type=int, default=100,
                    help="ignore HSPs shorter than this when counting")
    ap.add_argument("--clean-max", type=int, default=2,
                    help="tiles at or below this hit count are labelled clean")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("-o", "--out", type=Path, help="write TSV here instead of stdout")
    args = ap.parse_args(argv)

    if not shutil.which("blastn"):
        raise SystemExit("blastn not found on PATH")
    if args.tile < 50:
        raise SystemExit("--tile below 50 bp gives noise, not a map")

    _, seq = read_single_fasta(args.fasta)
    lo, hi = parse_region(args.region, len(seq))

    starts = list(range(lo, hi + 1, args.tile))
    with tempfile.TemporaryDirectory() as tmp:
        tiles = Path(tmp) / "tiles.fa"
        tiles.write_text("".join(
            f">t{s}\n{seq[s - 1:s - 1 + args.tile]}\n" for s in starts
        ))
        print(f"{len(starts)} tiles of {args.tile} bp over {lo}-{hi}", file=sys.stderr)
        proc = subprocess.run(
            ["blastn", "-task", "megablast", "-query", str(tiles), "-db", args.db,
             "-outfmt", "6 qseqid length", "-evalue", "1e-5",
             "-max_target_seqs", "5000", "-num_threads", str(args.threads)],
            capture_output=True, text=True,
        )
    if proc.returncode != 0:
        raise SystemExit(f"blastn exited {proc.returncode}: {proc.stderr.strip()[:300]}")

    counts: dict[int, int] = {s: 0 for s in starts}
    for line in proc.stdout.splitlines():
        q, ln = line.split("\t")
        if int(ln) >= args.min_hit_len:
            counts[int(q[1:])] += 1

    rows = ["start\tend\thits\tverdict"]
    for s in starts:
        n = counts[s]
        verdict = "clean" if n <= args.clean_max else ("mixed" if n <= 20 else "repeat")
        rows.append(f"{s}\t{min(s + args.tile - 1, hi)}\t{n}\t{verdict}")

    clean = sum(1 for s in starts if counts[s] <= args.clean_max)
    print(f"clean tiles: {clean}/{len(starts)}", file=sys.stderr)

    text = "\n".join(rows) + "\n"
    if args.out:
        args.out.write_text(text)
        print(f"wrote {args.out}", file=sys.stderr)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
