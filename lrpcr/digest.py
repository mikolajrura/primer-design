#!/usr/bin/env python3
"""Per-amplicon restriction digest notes, aware of a feature you must not cut.

Written for amplicons that carry something whose length is the readout: an
expanded trinucleotide tract, a variable insert. Two things matter then:

  * which enzymes cut *inside* that feature (never usable: they destroy the
    measurement), and
  * which enzymes release a small fragment containing the whole feature, so a
    40 kb product can still be sized on a short, well-resolved band.

Mark the protected span in the FASTA header as `protect=<start>-<end>`,
1-based and inclusive relative to the amplicon, or pass --protect for a single
amplicon. Without it the tool still reports cut sites and fragments.

Enzyme data comes from REBASE via Biopython's Bio.Restriction; the exact REBASE
release is whatever your Biopython bundles, and it is printed in each note.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

try:
    from Bio.Restriction import AllEnzymes, CommOnly, RestrictionBatch
    from Bio.Seq import Seq
except ImportError:  # pragma: no cover
    raise SystemExit("Biopython is required: pip install biopython")


def rebase_version() -> str:
    try:
        from Bio.Restriction import Restriction_Dictionary as RD
        src = Path(RD.__file__).read_text()[:2000]
        m = re.search(r"REBASE[^\n]*?(\d{3,4})[^\n]*?\((\d{4})\)", src)
        if m:
            return f"REBASE emboss v{m.group(1)} ({m.group(2)})"
    except Exception:
        pass
    return "REBASE release unknown"


def read_fasta(path: Path) -> list[tuple[str, str, str]]:
    """Return (name, header_rest, sequence) for every record."""
    out, name, rest, chunks = [], None, "", []
    for line in path.read_text().splitlines():
        line = line.rstrip()
        if line.startswith(">"):
            if name is not None:
                out.append((name, rest, "".join(chunks).upper()))
            parts = line[1:].split(None, 1)
            name, rest, chunks = parts[0], (parts[1] if len(parts) > 1 else ""), []
        elif name is not None:
            chunks.append(line.strip())
    if name is not None:
        out.append((name, rest, "".join(chunks).upper()))
    if not out:
        raise SystemExit(f"{path}: no FASTA records")
    return out


def protect_from_header(rest: str) -> tuple[int, int] | None:
    m = re.search(r"protect=(\d+)-(\d+)", rest)
    return (int(m.group(1)), int(m.group(2))) if m else None


def resolvable(frags: list[int]) -> bool:
    """Would every band be told apart on a normal agarose gel?"""
    s = sorted(frags, reverse=True)
    return all(a - b >= 500 or a / b >= 1.15 for a, b in zip(s, s[1:]))


def analyse(seq: str, protect, batch, min_fragment: int, max_cuts: int):
    length = len(seq)
    found = batch.search(Seq(seq))
    cutting = {str(e): sorted(v) for e, v in found.items() if v}
    nocut = sorted(str(e) for e, v in found.items() if not v)

    rows = []
    for enzyme, cuts in cutting.items():
        frags, prev = [], 0
        for c in cuts:
            frags.append(c - prev)
            prev = c
        frags.append(length - prev)

        in_prot = bool(protect) and any(protect[0] <= c <= protect[1] for c in cuts)
        frag_with = ""
        if protect and not in_prot:
            prev = 0
            for c in cuts + [length]:
                if prev < protect[0] and c >= protect[1]:
                    frag_with = c - prev
                    break
                prev = c
        diagnostic = (
            not in_prot
            and 1 <= len(cuts) <= max_cuts
            and min(frags) >= min_fragment
            and resolvable(frags)
        )
        rows.append({
            "enzyme": enzyme, "cuts": cuts, "frags": frags, "in_protected": in_prot,
            "frag_with_feature": frag_with, "min_frag": min(frags), "diagnostic": diagnostic,
        })
    rows.sort(key=lambda r: (not r["diagnostic"], len(r["cuts"]), -r["min_frag"]))
    return rows, nocut


def note(name: str, seq: str, protect, rows, nocut, n_enzymes: int,
         min_fragment: int, max_cuts: int) -> str:
    gc = 100 * sum(c in "GC" for c in seq) / len(seq)
    diag = [r for r in rows if r["diagnostic"]]
    once = [r for r in rows if len(r["cuts"]) == 1 and not r["in_protected"]]
    inprot = [r["enzyme"] for r in rows if r["in_protected"]]

    L = [f"DIGEST: {name}", "=" * 72,
         f"amplicon      : {len(seq)} bp, GC {gc:.2f}%",
         f"protected span: {f'{protect[0]}-{protect[1]}' if protect else 'none given'}",
         f"enzyme set    : {n_enzymes}, {rebase_version()}",
         f"cut / no cut  : {len(rows)} / {len(nocut)}", ""]

    L.append(f"[1] DIAGNOSTIC ({len(diag)}): outside the protected span, 1-{max_cuts} cuts,")
    L.append(f"    every fragment >= {min_fragment} bp and resolvable on a gel")
    if diag:
        L.append(f"    {'enzyme':<12}{'cuts':>5}  {'fragments (bp, desc)':<38}{'w/ feature':>11}")
        for r in diag[:25]:
            f = " + ".join(str(x) for x in sorted(r["frags"], reverse=True))
            if len(f) > 36:
                f = f[:33] + "..."
            L.append(f"    {r['enzyme']:<12}{len(r['cuts']):>5}  {f:<38}{r['frag_with_feature']:>11}")
        if len(diag) > 25:
            L.append(f"    ... {len(diag) - 25} more, see the TSV")
    else:
        L.append("    NONE")

    L += ["", f"[2] SINGLE CUTTERS outside the protected span ({len(once)}), for linearising"]
    L.append("    " + (", ".join(f"{r['enzyme']}({r['cuts'][0]})" for r in once[:30]) or "NONE"))
    if len(once) > 30:
        L.append(f"    ... {len(once) - 30} more")

    L += ["", f"[3] NON-CUTTERS ({len(nocut)}), free for cloning"]
    L.append("    " + ", ".join(nocut[:40]))
    if len(nocut) > 40:
        L.append(f"    ... {len(nocut) - 40} more, see the TSV")

    L += ["", f"[4] CUTTING INSIDE THE PROTECTED SPAN ({len(inprot)}): DO NOT USE for sizing"]
    L.append("    " + (", ".join(sorted(inprot)) if inprot else "none"))
    L += ["", "Cut positions are 1-based amplicon coordinates (top-strand cut is after"
              " the given base)."]
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Restriction digest notes per amplicon, protecting a named span.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--amplicons", type=Path, required=True,
                    help="FASTA file, or a directory of .fa files")
    ap.add_argument("--protect", help="span start-end applied to every amplicon; "
                                      "otherwise read from each header's protect=A-B")
    ap.add_argument("--min-fragment", type=int, default=250,
                    help="fragments below this are not reliably visible")
    ap.add_argument("--max-cuts", type=int, default=6, help="cuts still counted as diagnostic")
    ap.add_argument("--all-enzymes", action="store_true",
                    help="use every REBASE enzyme rather than commercially available only")
    ap.add_argument("-o", "--outdir", type=Path, required=True)
    args = ap.parse_args(argv)

    files = (sorted(args.amplicons.glob("*.fa")) if args.amplicons.is_dir()
             else [args.amplicons])
    if not files:
        raise SystemExit(f"{args.amplicons}: no .fa files")

    fixed = None
    if args.protect:
        a, b = args.protect.split("-")
        fixed = (int(a), int(b))

    enzymes = AllEnzymes if args.all_enzymes else CommOnly
    batch = RestrictionBatch(list(enzymes))
    args.outdir.mkdir(parents=True, exist_ok=True)

    combined = ["amplicon\tenzyme\tn_cuts\tcut_positions\tfragments_bp\t"
                "cuts_in_protected\tfragment_with_feature\tmin_fragment\tdiagnostic"]
    n_rec = 0
    for path in files:
        for name, rest, seq in read_fasta(path):
            protect = fixed or protect_from_header(rest)
            if protect and not (1 <= protect[0] < protect[1] <= len(seq)):
                raise SystemExit(f"{name}: protected span {protect} outside 1-{len(seq)}")
            rows, nocut = analyse(seq, protect, batch, args.min_fragment, args.max_cuts)
            (args.outdir / f"{name}__digest.txt").write_text(
                note(name, seq, protect, rows, nocut, len(enzymes),
                     args.min_fragment, args.max_cuts)
            )
            for r in rows:
                combined.append("\t".join(map(str, [
                    name, r["enzyme"], len(r["cuts"]),
                    ";".join(map(str, r["cuts"])),
                    ";".join(str(x) for x in sorted(r["frags"], reverse=True)),
                    "yes" if r["in_protected"] else "no",
                    r["frag_with_feature"], r["min_frag"],
                    "yes" if r["diagnostic"] else "",
                ])))
            n_rec += 1

    (args.outdir / "_digests.tsv").write_text("\n".join(combined) + "\n")
    print(f"{n_rec} amplicons -> {args.outdir} "
          f"({len(combined) - 1} rows in _digests.tsv)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
