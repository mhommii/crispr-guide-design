"""Step 4 — Search chromosome 17 for places each guide could also cut.

The problem
-----------
A guide is 20 bases long, but Cas9 tolerates mismatches. A guide designed
for TP53 may also bind a similar sequence elsewhere in the genome and cut
there. Those unintended cuts are off-target effects, and checking for them
is the reason guide design is not just "pick any 20-mer".

Tolerance is not uniform along the guide. Mismatches in the seed region -
roughly the 10-12 bases next to the PAM - usually prevent cutting.
Mismatches at the far (5') end are tolerated much more often. So a site
with 3 mismatches all at the PAM-distal end is a bigger worry than a site
with 2 mismatches inside the seed.

Why this file implements the search instead of calling Cas-OFFinder
------------------------------------------------------------------
Cas-OFFinder (Bae et al. 2014) is the standard tool and step 4 still writes
an input file for it (results/cas_offinder_input.txt). But it requires an
OpenCL runtime, and this machine has none registered, so it exits with
"No OpenCL devices found". Rather than skip the step, the same search is
done here directly. It is exhaustive over chromosome 17 - every NGG site on
both strands is compared against every guide - so it is not an
approximation of Cas-OFFinder, just a slower-to-write version of the same
idea without bulges.

How it is made fast enough
--------------------------
Each base is packed into 2 bits (A=00, C=01, G=10, T=11), so a 20-base
protospacer becomes one 40-bit integer. Comparing a guide to a site is then
an XOR: positions that differ leave a nonzero 2-bit pair. Folding each pair
down to a single bit and counting set bits gives the mismatch count for a
whole site in a few machine instructions, and numpy applies it to all ~10
million sites at once.

Limitation: chromosome 17 only, not the whole genome (the full reference is
~3 GB). A guide that looks clean here has been checked against about 2.7%
of the genome. Bulges (insertions/deletions between guide and DNA) are not
searched; Cas-OFFinder can do those.

Output:
  results/offtargets.tsv          every off-target site found
  results/offtarget_summary.tsv   per-guide counts by mismatch number
  results/cas_offinder_input.txt  input file to run the real tool elsewhere
"""

import gzip
import shutil
import tomllib
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
config = tomllib.loads((ROOT / "config.toml").read_text())
gene = config["gene"]

PROTOSPACER_LENGTH = 20
MAX_MISMATCHES = 4
SEED_LENGTH = 12          # bases next to the PAM where mismatches matter most
CHUNK = 4_000_000         # sites packed at a time, to keep memory modest

data, results = ROOT / "data", ROOT / "results"

# --- Load the reference ----------------------------------------------------
genome_dir = data / "genome"
genome_dir.mkdir(exist_ok=True)
chromosome_path = genome_dir / "chr17.fa"
if not chromosome_path.exists():
    archive = data / "chr17.fa.gz"
    if not archive.exists():
        raise SystemExit("data/chr17.fa.gz missing - see README for the download.")
    print("Decompressing chr17.fa.gz ...")
    with gzip.open(archive, "rb") as source, chromosome_path.open("wb") as target:
        shutil.copyfileobj(source, target)

print("Loading chr17 ...")
with chromosome_path.open("rb") as handle:
    handle.readline()  # drop the ">chr17" header
    sequence = handle.read().replace(b"\n", b"").upper()

# Encode bases as 0-3; anything else (N, masked repeats already uppercased)
# becomes 255 so we can exclude those windows.
lookup = np.full(256, 255, dtype=np.uint8)
for base, code in zip(b"ACGT", range(4)):
    lookup[base] = code
codes = lookup[np.frombuffer(sequence, dtype=np.uint8)]
del sequence
print(f"chr17: {len(codes):,} bp  ({np.count_nonzero(codes == 255):,} are N)")

WEIGHTS = (4 ** np.arange(PROTOSPACER_LENGTH - 1, -1, -1)).astype(np.uint64)
SEED_MASK = np.uint64((1 << (2 * SEED_LENGTH)) - 1)
PAIR_MASK = np.uint64(0x5555555555555555)


def pack_windows(starts: np.ndarray, source: np.ndarray) -> np.ndarray:
    """Pack the 20-base window beginning at each start index into a uint64."""
    packed = np.empty(len(starts), dtype=np.uint64)
    for begin in range(0, len(starts), CHUNK):
        block = starts[begin:begin + CHUNK]
        window = source[block[:, None] + np.arange(PROTOSPACER_LENGTH)]
        packed[begin:begin + CHUNK] = window.astype(np.uint64) @ WEIGHTS
    return packed


def find_sites(source: np.ndarray, strand: str):
    """Every NGG site on one strand: packed protospacer + start coordinate."""
    # PAM is NGG at positions p, p+1, p+2 -> protospacer occupies p-20 .. p-1.
    # Element k of (is_g[1:-1] & is_g[2:]) is True when source[k+1] and
    # source[k+2] are both G, which means the PAM starts at k itself.
    is_g = source == 2
    pam_starts = np.flatnonzero(is_g[1:-1] & is_g[2:])
    pam_starts = pam_starts[pam_starts >= PROTOSPACER_LENGTH]
    protospacer_starts = pam_starts - PROTOSPACER_LENGTH

    # Drop any window (protospacer + PAM) containing an ambiguous base.
    valid = np.ones(len(protospacer_starts), dtype=bool)
    for begin in range(0, len(protospacer_starts), CHUNK):
        block = protospacer_starts[begin:begin + CHUNK]
        window = source[block[:, None] + np.arange(PROTOSPACER_LENGTH + 3)]
        valid[begin:begin + CHUNK] = ~(window == 255).any(axis=1)
    protospacer_starts = protospacer_starts[valid]

    print(f"  {strand} strand: {len(protospacer_starts):,} NGG sites")
    return pack_windows(protospacer_starts, source), protospacer_starts


print("Indexing NGG sites on both strands ...")
forward_packed, forward_starts = find_sites(codes, "+")

# Reverse complement: for codes A=0,C=1,G=2,T=3 the complement is 3 - code.
reverse_codes = 3 - codes[::-1]
reverse_codes[codes[::-1] == 255] = 255
reverse_packed, reverse_starts = find_sites(reverse_codes, "-")
chromosome_length = len(codes)
del codes

site_packed = np.concatenate([forward_packed, reverse_packed])
site_start = np.concatenate([forward_starts, reverse_starts])
site_strand = np.concatenate([
    np.zeros(len(forward_packed), dtype=np.uint8),
    np.ones(len(reverse_packed), dtype=np.uint8),
])
del forward_packed, reverse_packed, forward_starts, reverse_starts
print(f"  total: {len(site_packed):,} sites indexed")


def encode(protospacer: str) -> np.uint64:
    value = 0
    for base in protospacer:
        value = (value << 2) | "ACGT".index(base)
    return np.uint64(value)


def decode(value: int, length: int = PROTOSPACER_LENGTH) -> str:
    return "".join("ACGT"[(value >> (2 * shift)) & 0b11]
                   for shift in range(length - 1, -1, -1))


# --- Search ----------------------------------------------------------------
guides = pd.read_csv(results / "filtered_guides.tsv", sep="\t")
print(f"\nSearching {len(guides)} guides against chr17 "
      f"(up to {MAX_MISMATCHES} mismatches) ...")

records = []
for count, protospacer in enumerate(guides["protospacer"], start=1):
    difference = site_packed ^ encode(protospacer)
    # Collapse each 2-bit pair to one bit: set if the two bases differ.
    folded = (difference | (difference >> np.uint64(1))) & PAIR_MASK
    mismatches = np.bitwise_count(folded)

    hit_index = np.flatnonzero(mismatches <= MAX_MISMATCHES)
    seed_mismatches = np.bitwise_count(folded[hit_index] & SEED_MASK)

    for position, seed in zip(hit_index, seed_mismatches):
        strand = "-" if site_strand[position] else "+"
        start = int(site_start[position])
        # Map reverse-strand coordinates back onto the forward chromosome.
        if strand == "-":
            start = chromosome_length - start - (PROTOSPACER_LENGTH + 3)
        records.append({
            "protospacer": protospacer,
            "site": decode(int(site_packed[position])),
            "chromosome": "chr17",
            "start": start,
            "strand": strand,
            "mismatches": int(mismatches[position]),
            "seed_mismatches": int(seed),
        })
    if count % 50 == 0:
        print(f"  {count}/{len(guides)} guides searched")

hits = pd.DataFrame(records)
hits.sort_values(["protospacer", "mismatches"], inplace=True)
hits.to_csv(results / "offtargets.tsv", sep="\t", index=False)
print(f"\n{len(hits):,} sites found in total (including each guide's own target).")

# Self-check: every guide was taken from TP53, which is on chr17, so every
# guide must match itself somewhere in this chromosome. If any guide has no
# perfect match the indexing is wrong, and the off-target counts below would
# be meaningless.
perfect = hits[hits["mismatches"] == 0]
without_perfect_match = set(guides["protospacer"]) - set(perfect["protospacer"])
if without_perfect_match:
    raise SystemExit(
        f"Indexing error: {len(without_perfect_match)} guides have no perfect "
        f"match on chr17, but every guide came from TP53 on chr17.\n"
        f"Example: {sorted(without_perfect_match)[0]}"
    )
print(f"Self-check passed: all {len(guides)} guides found their own target site.")

# --- Per-guide summary -----------------------------------------------------
rows = []
for protospacer in guides["protospacer"]:
    guide_hits = hits[hits["protospacer"] == protospacer]
    counts = guide_hits["mismatches"].value_counts()
    off_targets = guide_hits[guide_hits["mismatches"] > 0]
    rows.append({
        "protospacer": protospacer,
        "perfect_matches": int(counts.get(0, 0)),
        "mm1": int(counts.get(1, 0)),
        "mm2": int(counts.get(2, 0)),
        "mm3": int(counts.get(3, 0)),
        "mm4": int(counts.get(4, 0)),
        "offtargets_total": int(len(off_targets)),
        # The worrying ones: few mismatches, and none of them in the seed.
        "offtargets_intact_seed": int((off_targets["seed_mismatches"] == 0).sum()),
    })

summary = pd.DataFrame(rows)
summary.to_csv(results / "offtarget_summary.tsv", sep="\t", index=False)

# --- Also write a Cas-OFFinder input file for use on an OpenCL machine -----
with (results / "cas_offinder_input.txt").open("w", encoding="utf-8") as handle:
    # Relative to the repository root, so the file is portable.
    handle.write(f"{genome_dir.relative_to(ROOT).as_posix()}\n"
                 f"{'N' * PROTOSPACER_LENGTH}NGG\n")
    for protospacer in guides["protospacer"]:
        handle.write(f"{protospacer}NNN {MAX_MISMATCHES}\n")

print(f"Guides whose only perfect chr17 match is one site: "
      f"{(summary['perfect_matches'] == 1).sum()} of {len(summary)}")
print(f"Guides with no off-target at 1-2 mismatches:       "
      f"{((summary['mm1'] + summary['mm2']) == 0).sum()}")
print(f"Median off-targets per guide (1-4 mismatches):     "
      f"{summary['offtargets_total'].median():.0f}")
print("\nSaved: results/offtargets.tsv, results/offtarget_summary.tsv, "
      "results/cas_offinder_input.txt")
