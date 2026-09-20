# Method

The order matters as much as the individual checks. Running them in the wrong
order is how an afternoon disappears.

## Order of operations

```
1. repeat map of the target region        cheap, kills whole windows at once
2. primer3 on the surviving windows       generates candidates
3. local filters (3' end, complexity)     free, no BLAST
4. BLAST vs the intended template         must be 100% / 0 mm / 0 gaps
5. 3'-anchor screen                       the expensive one, on a shortlist
6. digest planning, if a tract is sized
```

Steps 3 and 4 are cheap and cut the candidate set by an order of magnitude.
Do them **before** step 5.

## Why step 1 first

Designing first and filtering afterwards means running primer3 in windows
where nothing can succeed. In one 31 kb region, half the tiles were inside
repeat families; every candidate from those tiles had tens of thousands of
full-length genomic hits. A tiling pass costs one megablast run and tells you
where to point primer3.

## Why 5 is last, and chunked

BLASTing short exact *k*-mers against a mammalian genome is memory-hungry: a
handful of 12-mer queries can return hundreds of thousands of HSPs, and the
result set is held in memory. Screening ~1000 candidates at once will exhaust
RAM. `anchors.py` chunks queries (`--chunk`, default 8) and retries a failing
chunk one query at a time, but the real fix is to shortlist first.

If you must run it on a big set, cap memory so a failure is a failure and not
a dead machine:

```bash
( ulimit -v 6000000; python3 lrpcr/anchors.py ... )
```

## Local filters worth applying before any BLAST

Two cost nothing and remove a lot:

**3'-end composition.** No more than 3 G/C in the last 5 nt. A GC-clamped
3' end is standard advice for stability, but a *run* of G/C there is also what
makes a promiscuous anchor.

**Sequence complexity.** Reject candidates whose distinct-3-mer fraction is
below ~0.9. This catches microsatellites: a `(GT)n` primer can be unique at
full length and still have five figures of partial hits.

```python
def complexity(s):
    k = [s[i:i+3] for i in range(len(s) - 2)]
    return len(set(k)) / len(k)
```

## Interpreting the anchor numbers

There is no universal threshold; the numbers are for ranking within one design
round on one genome. What held up in practice:

- `self_products` tracks laddering better than `anchors` does. A primer with
  many anchors but no opposing pairs within range has nowhere to make a
  product on its own.
- Compare against primers you have already run. If a primer that laddered on a
  gel scores 8-20 self-products and a new candidate scores 0-2, that is a
  meaningful ranking. An absolute cutoff is not.
- Full-length uniqueness and anchor count are independent. Check both.

**The control that settles it:** run the suspect primer alone, without a
partner. If bands appear, the self-product model was right. Nothing in this
repository substitutes for that lane.

## Digest planning for a sized tract

When the product is long, the length difference you care about is a small
fraction of it and will not resolve. Look for an enzyme that releases a short
fragment containing the whole tract; then the measurement moves onto a band
where it is legible, and the amplicon length stops mattering.

What `digest.py` reports per amplicon:

- **diagnostic**: outside the protected span, few cuts, every fragment above
  `--min-fragment` and mutually resolvable (at least 500 bp apart, or 15% different)
- **single cutters** outside the span, for linearising
- **non-cutters**, free for cloning
- **cutting inside the protected span**: these destroy the readout silently,
  which is the whole reason the section exists

Two things the tool does not know:

- **Methylation.** Cut maps hold on PCR product. On genomic DNA, a
  CpG-sensitive enzyme in a GC-rich region may not cut.
- **Availability.** REBASE lists enzymes as commercially available that may in
  practice come from one niche supplier. Check before building a protocol on
  a rare enzyme, however clean its map looks.

## Reproducing the numbers

Every count in the README came from the demo script against GRCm39. Nothing is
estimated. If a number here cannot be reproduced by a command, it is a bug in
the documentation.
