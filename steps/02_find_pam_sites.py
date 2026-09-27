"""Step 2 — Find every SpCas9 target site in the gene.

How SpCas9 chooses where to cut
-------------------------------
Cas9 only cuts where it finds a short DNA motif called a PAM
(protospacer adjacent motif). For the common Streptococcus pyogenes Cas9
the PAM is 5'-NGG-3': any base, then two Gs.

The 20 bases immediately 5' of the PAM are the protospacer. The guide RNA
carries that same 20-base sequence, so it pairs with the opposite strand
and positions Cas9 to cut between positions 17 and 18 of the protospacer
(3 bp upstream of the PAM).

    5'-...NNNNNNNNNNNNNNNNNNNN NGG...-3'
           protospacer (20 nt)  PAM
                          ^ blunt cut, 3 bp before the PAM

DNA is double-stranded, and Cas9 does not care which strand carries the
PAM. So we search the forward strand for NGG and the forward strand for
CCN, which is an NGG read on the reverse strand. Every hit on the reverse
strand gives a protospacer that is the reverse complement of the forward
sequence.

Output:
  results/all_sites.tsv  every candidate site with coordinates and strand
"""

import re
import tomllib
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parent.parent
config = tomllib.loads((ROOT / "config.toml").read_text())
gene = config["gene"]

record = SeqIO.read(ROOT / "data" / f"{gene}.gb", "genbank")
sequence = str(record.seq).upper()

PROTOSPACER_LENGTH = 20

# Coding exons, so we can tell which guides land in protein-coding sequence.
# A gene has several annotated transcripts; we take the union of their exons.
coding_ranges = []
for feature in record.features:
    if feature.type == "CDS":
        for part in feature.location.parts:
            coding_ranges.append((int(part.start), int(part.end)))


def in_coding_exon(start: int, end: int) -> bool:
    """True if the site overlaps any annotated coding exon."""
    return any(start < exon_end and end > exon_start
               for exon_start, exon_end in coding_ranges)


sites = []

# --- Forward strand: look for NGG ------------------------------------------
# The regex uses a lookahead so that overlapping PAMs are all reported
# (GGG contains two valid NGG motifs).
for match in re.finditer(r"(?=([ACGT]GG))", sequence):
    pam_start = match.start()
    protospacer_start = pam_start - PROTOSPACER_LENGTH
    if protospacer_start < 0:
        continue
    protospacer = sequence[protospacer_start:pam_start]
    if "N" in protospacer:
        continue
    sites.append({
        "protospacer": protospacer,
        "pam": match.group(1),
        "strand": "+",
        "start": protospacer_start,
        "end": pam_start,
        # Cas9 cuts 3 bp 5' of the PAM
        "cut_site": pam_start - 3,
        "in_cds": in_coding_exon(protospacer_start, pam_start),
    })

# --- Reverse strand: an NGG there reads as CCN here -------------------------
for match in re.finditer(r"(?=(CC[ACGT]))", sequence):
    pam_start = match.start()
    protospacer_end = pam_start + 3 + PROTOSPACER_LENGTH
    if protospacer_end > len(sequence):
        continue
    forward_slice = sequence[pam_start + 3:protospacer_end]
    if "N" in forward_slice:
        continue
    sites.append({
        # The guide matches the reverse strand, so reverse complement it.
        "protospacer": str(Seq(forward_slice).reverse_complement()),
        "pam": str(Seq(match.group(1)).reverse_complement()),
        "strand": "-",
        "start": pam_start + 3,
        "end": protospacer_end,
        "cut_site": pam_start + 3 + 3,
        "in_cds": in_coding_exon(pam_start + 3, protospacer_end),
    })

sites.sort(key=lambda s: (s["start"], s["strand"]))

results_dir = ROOT / "results"
results_dir.mkdir(exist_ok=True)
output = results_dir / "all_sites.tsv"

columns = ["protospacer", "pam", "strand", "start", "end", "cut_site", "in_cds"]
with output.open("w", encoding="utf-8") as handle:
    handle.write("\t".join(columns) + "\n")
    for site in sites:
        handle.write("\t".join(str(site[column]) for column in columns) + "\n")

forward = sum(1 for s in sites if s["strand"] == "+")
coding = sum(1 for s in sites if s["in_cds"])
print(f"Sequence searched: {len(sequence):,} bp of {gene}")
print(f"Coding exon blocks: {len(coding_ranges)} (union of all transcripts)")
print(f"Candidate sites:   {len(sites):,}  ({forward:,} on + strand, "
      f"{len(sites) - forward:,} on - strand)")
print(f"In a coding exon:  {coding:,}")
print(f"Density:           {len(sites) / len(sequence):.2f} sites per bp")
print(f"Saved:             results/all_sites.tsv")
