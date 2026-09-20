#!/usr/bin/env python3
"""Offline checks. No genome or BLAST database needed.

Run:  python3 tests/test_smoke.py
"""
from __future__ import annotations

import random
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from lrpcr import anchors, digest  # noqa: E402

FAILS: list[str] = []


def check(label: str, got, want):
    if got == want:
        print(f"  ok   {label}")
    else:
        print(f"  FAIL {label}: got {got!r}, want {want!r}")
        FAILS.append(label)


def synthetic_amplicon(n_repeats: int = 40, seed: int = 7) -> tuple[str, int, int]:
    """Random flanks around a (CTG)n tract; returns seq and the tract span."""
    rng = random.Random(seed)
    flank = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    left, right = flank(600), flank(900)
    tract = "CTG" * n_repeats
    seq = left + tract + right
    return seq, len(left) + 1, len(left) + len(tract)


def test_fasta_roundtrip():
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "p.fa"
        p.write_text(">a desc here\nACGT\nacgt\n>b\nTTTT\n")
        d = anchors.read_fasta(p)
        check("read_fasta names", sorted(d), ["a", "b"])
        check("read_fasta uppercases and joins", d["a"], "ACGTACGT")


def test_products_geometry():
    """One plus anchor at 100 and one minus anchor at 1100 make one product."""
    idx = anchors.index([
        ("F", "chr1", 100, "plus"),
        ("R", "chr1", 1100, "minus"),
        ("R", "chr1", 90, "minus"),      # upstream of F: not a product
        ("R", "chr2", 500, "minus"),     # other contig: not a product
    ])
    check("product within range", anchors.products(idx, "F", "R", 5000), 1)
    check("product beyond max distance", anchors.products(idx, "F", "R", 500), 0)
    check("unknown primer", anchors.products(idx, "F", "ZZZ", 5000), 0)


def test_resolvable():
    check("clearly separated", digest.resolvable([5000, 1000]), True)
    check("near-identical bands", digest.resolvable([1000, 990]), False)
    check("far apart in bp", digest.resolvable([9000, 8000]), True)


def test_digest_protects_tract():
    seq, lo, hi = synthetic_amplicon()
    from Bio.Restriction import CommOnly, RestrictionBatch
    batch = RestrictionBatch(list(CommOnly))
    rows, nocut = digest.analyse(seq, (lo, hi), batch, min_fragment=250, max_cuts=6)
    check("every enzyme is either cutting or not", len(rows) + len(nocut), len(CommOnly))

    for r in rows:
        total = sum(r["frags"])
        if total != len(seq):
            check(f"fragments of {r['enzyme']} sum to amplicon length", total, len(seq))
            break
    else:
        print("  ok   fragment lengths always sum to the amplicon")

    bad = [r["enzyme"] for r in rows if r["in_protected"] and r["diagnostic"]]
    check("no enzyme cutting the tract is called diagnostic", bad, [])

    for r in rows:
        if r["frag_with_feature"]:
            if not (r["frag_with_feature"] >= hi - lo + 1):
                check(f"{r['enzyme']} feature fragment spans the tract", False, True)
                break
    else:
        print("  ok   reported feature fragments cover the tract")


def test_cli_end_to_end():
    seq, lo, hi = synthetic_amplicon()
    with tempfile.TemporaryDirectory() as tmp:
        fa = Path(tmp) / "amp.fa"
        fa.write_text(f">demo protect={lo}-{hi}\n" +
                      "\n".join(seq[i:i + 70] for i in range(0, len(seq), 70)) + "\n")
        out = Path(tmp) / "out"
        rc = subprocess.run(
            [sys.executable, str(ROOT / "lrpcr" / "digest.py"),
             "--amplicons", str(fa), "-o", str(out)],
            capture_output=True, text=True,
        )
        check("digest.py exit code", rc.returncode, 0)
        check("note written", (out / "demo__digest.txt").exists(), True)
        check("combined tsv written", (out / "_digests.tsv").exists(), True)
        text = (out / "demo__digest.txt").read_text()
        check("note states the protected span", f"{lo}-{hi}" in text, True)


if __name__ == "__main__":
    for fn in (test_fasta_roundtrip, test_products_geometry, test_resolvable,
               test_digest_protects_tract, test_cli_end_to_end):
        print(fn.__name__)
        fn()
    print()
    if FAILS:
        print(f"{len(FAILS)} FAILED: {', '.join(FAILS)}")
        raise SystemExit(1)
    print("all checks passed")
