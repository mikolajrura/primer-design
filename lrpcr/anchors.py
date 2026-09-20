#!/usr/bin/env python3
"""3'-anchor mispriming screen for PCR primers.

Full-length BLAST uniqueness is a weak predictor of clean amplification. What
actually seeds a spurious product is the primer's 3' end: a short perfect match
there is extendable even when the rest of the primer does not pair. This module
counts those 3' anchors genome-wide and, more usefully, counts *pairs* of
opposing anchors close enough to yield a product.

Two numbers come out per primer:

  anchors         how many places the 3'-terminal k-mer matches perfectly
  self_products   pairs of opposing anchors of the SAME primer within
                  --max-product, i.e. products one primer can make unaided

A primer that is unique at full length can still have a thousand anchors and a
dozen self-products. Those are the ones that ladder on a gel.

BLAST of short exact k-mers against a mammalian genome is memory-hungry, so
queries are chunked and each chunk runs as its own subprocess; see --chunk.
"""
from __future__ import annotations

import argparse
import bisect
import collections
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

Hit = tuple[str, str, int, str]  # query, subject, leftmost pos, strand


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        elif name is not None:
            seqs[name] += line.upper()
    if not seqs:
        raise SystemExit(f"{path}: no FASTA records")
    return seqs


def blast_anchors(fasta: Path, db: str, k: int, threads: int) -> list[Hit]:
    """Every exact k-mer hit of the queries in `fasta` against `db`."""
    cmd = [
        "blastn", "-task", "blastn", "-word_size", str(k), "-perc_identity", "100",
        "-query", str(fasta), "-db", db,
        "-outfmt", "6 qseqid sseqid sstart send length",
        "-evalue", "10000", "-max_target_seqs", "500000",
        "-num_threads", str(threads), "-dust", "no",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"blastn exited {proc.returncode}: {proc.stderr.strip()[:300]}")
    hits: list[Hit] = []
    for line in proc.stdout.splitlines():
        q, s, a, b, ln = line.split("\t")
        if int(ln) != k:  # BLAST can report shorter HSPs even at word_size k
            continue
        a, b = int(a), int(b)
        hits.append((q, s, min(a, b), "plus" if a < b else "minus"))
    return hits


def collect(primers: dict[str, str], db: str, k: int, chunk: int, threads: int) -> list[Hit]:
    names = list(primers)
    hits: list[Hit] = []
    with tempfile.TemporaryDirectory() as tmp:
        for i in range(0, len(names), chunk):
            part = names[i:i + chunk]
            fa = Path(tmp) / f"chunk{i}.fa"
            fa.write_text("".join(f">{n}\n{primers[n][-k:]}\n" for n in part))
            try:
                hits += blast_anchors(fa, db, k, threads)
            except RuntimeError as exc:
                # One oversized query should not cost the whole chunk.
                print(f"  chunk at {i}: {exc}\n  retrying one at a time", file=sys.stderr)
                for n in part:
                    one = Path(tmp) / "one.fa"
                    one.write_text(f">{n}\n{primers[n][-k:]}\n")
                    try:
                        hits += blast_anchors(one, db, k, 1)
                    except RuntimeError as exc2:
                        print(f"  {n}: FAILED, counted as 0 anchors. {exc2}", file=sys.stderr)
    return hits


def index(hits: list[Hit]):
    idx = collections.defaultdict(
        lambda: collections.defaultdict(lambda: {"plus": [], "minus": []})
    )
    for q, s, pos, strand in hits:
        idx[q][s][strand].append(pos)
    for q in idx:
        for s in idx[q]:
            for st in idx[q][s]:
                idx[q][s][st].sort()
    return idx


def products(idx, fwd: str, rev: str, max_product: int) -> int:
    """Anchor pairs that could prime a product shorter than `max_product`.

    An anchor whose k-mer matches the plus strand acts as a forward primer
    there; its partner must match the minus strand downstream, and vice versa.
    """
    if fwd not in idx or rev not in idx:
        return 0
    total = 0
    for contig, by_strand in idx[fwd].items():
        if contig not in idx[rev]:
            continue
        for f_strand, r_strand in (("plus", "minus"), ("minus", "plus")):
            starts, ends = by_strand[f_strand], idx[rev][contig][r_strand]
            if not starts or not ends:
                continue
            for p in starts:
                if f_strand == "plus":
                    total += bisect.bisect_left(ends, p + max_product) - bisect.bisect_right(ends, p)
                else:
                    total += bisect.bisect_left(ends, p) - bisect.bisect_left(ends, p - max_product)
    return total


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Count 3' anchors and potential spurious products for PCR primers.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--primers", type=Path, required=True, help="FASTA of primer sequences")
    ap.add_argument("--db", required=True, help="BLAST nucleotide database (makeblastdb prefix)")
    ap.add_argument("--anchor-len", type=int, default=12,
                    help="length of the 3'-terminal k-mer treated as the anchor")
    ap.add_argument("--max-product", type=int, default=60000,
                    help="longest product the polymerase could plausibly make")
    ap.add_argument("--chunk", type=int, default=8,
                    help="primers per blastn call; lower it if BLAST runs out of memory")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--pairs", type=Path,
                    help="optional TSV of forward<TAB>reverse pairs to score as well")
    ap.add_argument("-o", "--out", type=Path, help="write TSV here instead of stdout")
    args = ap.parse_args(argv)

    if not shutil.which("blastn"):
        raise SystemExit("blastn not found on PATH")
    if args.anchor_len < 8:
        raise SystemExit("--anchor-len below 8 gives hit counts dominated by chance")

    primers = read_fasta(args.primers)
    short = [n for n, s in primers.items() if len(s) < args.anchor_len]
    if short:
        raise SystemExit(f"shorter than --anchor-len: {', '.join(short)}")

    print(f"{len(primers)} primers, anchor = last {args.anchor_len} nt", file=sys.stderr)
    idx = index(collect(primers, args.db, args.anchor_len, args.chunk, args.threads))

    rows = ["primer\tanchors\tself_products"]
    for n in primers:
        n_anchor = sum(len(v[st]) for v in idx.get(n, {}).values() for st in ("plus", "minus"))
        rows.append(f"{n}\t{n_anchor}\t{products(idx, n, n, args.max_product)}")

    if args.pairs:
        rows += ["", "forward\treverse\tpair_products"]
        for line in args.pairs.read_text().splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            f, r = line.split("\t")[:2]
            rows.append(f"{f}\t{r}\t{products(idx, f, r, args.max_product)}")

    text = "\n".join(rows) + "\n"
    if args.out:
        args.out.write_text(text)
        print(f"wrote {args.out}", file=sys.stderr)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
