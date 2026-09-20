# primer-design

Three small tools for designing and vetting **long-range PCR primers across a
repeat expansion**, where the length of the product is the measurement.

They exist because the usual checks were not enough. Primers that were unique
at full length against the genome, had matched Tm, no hairpin and no
heterodimer still produced several bands on a gel. The three screens here look
at what those checks miss.

```
lrpcr/anchors.py     3' anchors and self-products    why a "unique" primer ladders
lrpcr/repeatmap.py   repeat density before design    where not to bother designing
lrpcr/digest.py      digest notes per amplicon       sizing a tract in a 40 kb product
```

Nothing here replaces primer3. It runs after it, on the candidates primer3
returns, and before you spend money on synthesis.

## Install

```bash
python3 -m venv .venv && ./.venv/bin/pip install -r requirements.txt
```

External tools, expected on `PATH`: `blastn` and `makeblastdb` (NCBI BLAST+),
plus `samtools` for the demo. A BLAST database of your genome:

```bash
makeblastdb -in GRCm39.fa -dbtype nucl -out ~/genomes/blastdb/GRCm39
```

Offline checks, no genome needed:

```bash
./.venv/bin/python tests/test_smoke.py
```

## 1. `anchors.py`: why a unique primer still ladders

Full-length uniqueness is the wrong question. What seeds a spurious product is
the primer's **3' end**: a short perfect match there extends even when the rest
of the primer does not pair. The tool BLASTs the 3'-terminal *k*-mer (12 nt by
default), then counts pairs of opposing anchors close enough to make a product.

`self_products` is the number a **single primer, with no partner**, could make.

```bash
python3 lrpcr/anchors.py --primers primers.fa --db ~/genomes/blastdb/GRCm39 \
        --anchor-len 12 --max-product 60000 -o anchors.tsv
```

Add `--pairs pairs.tsv` (two columns, forward and reverse) to score real pairs
as well.

There is a decisive bench control for what this predicts: **run a suspect
primer alone, without its partner.** Bands means the model was right.

## 2. `repeatmap.py`: look before you design

Designing first and discarding repetitive candidates afterwards wastes most of
the work: kilobase-scale windows of a mammalian locus sit inside SINE/LINE
families, and *every* candidate from such a window fails. Tile the region
first.

```bash
python3 lrpcr/repeatmap.py --fasta locus.fa --db ~/genomes/blastdb/GRCm39 \
        --region 99000-130000 --tile 500 -o repeat_map.tsv
```

A unique tile scores 1 (itself). A tile inside a B1 element scores tens of
thousands. Design only in the clean ones.

## 3. `digest.py`: measuring a tract inside a long product

A 798 bp difference is invisible on a 40 kb band. But if a digest releases a
short fragment containing the whole tract, the **small band** carries the
measurement and the product length stops mattering.

Mark the span that must survive in the FASTA header:

```
>amp_1 protect=68-865
```

```bash
python3 lrpcr/digest.py --amplicons amplicons/ --min-fragment 250 -o digests/
```

Each amplicon gets a note with four sections: diagnostic enzymes, single
cutters, non-cutters, and **enzymes that cut inside the protected span**,
the ones that silently destroy the readout. Plus one combined TSV.

Enzymes come from REBASE via Biopython; the release is printed in every note.

## Demo, on public reference only

```bash
bash examples/demo_grcm39.sh
```

Pulls 30 kb of GRCm39, maps repeat density, and screens three 30-mers taken
straight off the reference. Observed on GRCm39 (mouse, C57BL/6J):

| primer | anchors | self_products |
|---|---|---|
| `demo_at_400` | 496 | 2 |
| `demo_at_2400` | 40 | 0 |
| `demo_at_3400` | 3662 | 350 |

with 15 of 31 tiles clean, 15 repeat, 1 mixed.

Note that `demo_at_2400` comes from a tile the map calls *repeat* yet has the
fewest anchors, while `demo_at_3400` sits in a *clean* tile and has the most.
**The two signals are independent**, so a candidate has to pass both.

## Caveats

- Everything here is computed, not measured. It narrows what you order; it
  does not tell you a PCR will work.
- The anchor model is a heuristic: a perfect *k*-mer at the 3' end, nothing
  about mismatch tolerance, secondary structure or polymerase behaviour. It
  ranks candidates, it does not predict products.
- Off-target counts are only as good as the assembly you BLAST against. On a
  non-reference strain background they are indicative, not verified.
- `digest.py` ignores methylation. Cut maps hold for PCR product; on genomic
  DNA a CpG-sensitive enzyme may not cut at all.

## Licence

MIT.
