"""Known-answer test for the packed mismatch counting used in step 4.

Step 4 compares a guide to millions of sites at once by packing each
20-base sequence into a 40-bit integer and counting differing 2-bit pairs.
That is fast but easy to get subtly wrong, so this checks it against a
plain, obviously-correct Python comparison.

Run:  .venv\\Scripts\\python.exe tests\\test_mismatch_counting.py
"""

import random
import sys

import numpy as np

PROTOSPACER_LENGTH = 20
SEED_LENGTH = 12
SEED_MASK = np.uint64((1 << (2 * SEED_LENGTH)) - 1)
PAIR_MASK = np.uint64(0x5555555555555555)


def encode(protospacer: str) -> np.uint64:
    value = 0
    for base in protospacer:
        value = (value << 2) | "ACGT".index(base)
    return np.uint64(value)


def packed_mismatches(guide: str, site: str):
    """The fast path used in step 4."""
    difference = encode(site) ^ encode(guide)
    folded = (difference | (difference >> np.uint64(1))) & PAIR_MASK
    return int(np.bitwise_count(folded)), int(np.bitwise_count(folded & SEED_MASK))


def plain_mismatches(guide: str, site: str):
    """The obvious slow version: compare base by base."""
    total = sum(1 for a, b in zip(guide, site) if a != b)
    seed = sum(1 for a, b in zip(guide[-SEED_LENGTH:], site[-SEED_LENGTH:]) if a != b)
    return total, seed


failures = []


def check(description, got, expected):
    if got == expected:
        print(f"  PASS  {description}")
    else:
        print(f"  FAIL  {description}: expected {expected!r}, got {got!r}")
        failures.append(description)


print("Test 1: identical sequences")
guide = "ACGTACGTACGTACGTACGT"
check("no mismatches against itself", packed_mismatches(guide, guide), (0, 0))

print("\nTest 2: a single mismatch at the PAM-distal end")
# Position 0 is 20 bases from the PAM, so it is outside the 12-nt seed.
distal = "C" + guide[1:]
check("counted as 1 mismatch, 0 in seed", packed_mismatches(guide, distal), (1, 0))

print("\nTest 3: a single mismatch inside the seed")
# Last base sits right against the PAM, well inside the seed.
seed_variant = guide[:-1] + ("A" if guide[-1] != "A" else "C")
check("counted as 1 mismatch, 1 in seed",
      packed_mismatches(guide, seed_variant), (1, 1))

print("\nTest 4: completely different sequences")
check("all 20 positions differ",
      packed_mismatches("A" * 20, "T" * 20), (20, 12))

print("\nTest 5: agrees with base-by-base counting on 20,000 random pairs")
random.seed(0)
mismatch_disagreements = 0
seed_disagreements = 0
for _ in range(20_000):
    a = "".join(random.choice("ACGT") for _ in range(PROTOSPACER_LENGTH))
    b = "".join(random.choice("ACGT") for _ in range(PROTOSPACER_LENGTH))
    fast_total, fast_seed = packed_mismatches(a, b)
    slow_total, slow_seed = plain_mismatches(a, b)
    mismatch_disagreements += fast_total != slow_total
    seed_disagreements += fast_seed != slow_seed
check("mismatch counts always agree", mismatch_disagreements, 0)
check("seed counts always agree", seed_disagreements, 0)

print("\nTest 6: reverse complement encoding (3 - code)")
# Step 4 builds the reverse strand with `3 - code`. Check that really is
# the complement for the A=0, C=1, G=2, T=3 scheme.
codes = {"A": 0, "C": 1, "G": 2, "T": 3}
complements = {"A": "T", "C": "G", "G": "C", "T": "A"}
check("3 - code is the complement for every base",
      all(3 - codes[base] == codes[complements[base]] for base in codes), True)

print("\nTest 7: PAM positions found in a packed genome array")
# This is the indexing used in step 4 to locate NGG sites on a chromosome.
# An off-by-one here silently shifts every protospacer by one base, so it
# is checked against a sequence where the answer can be counted by eye.


def find_pam_starts(codes: np.ndarray) -> np.ndarray:
    is_g = codes == 2
    return np.flatnonzero(is_g[1:-1] & is_g[2:])


def to_codes(sequence: str) -> np.ndarray:
    return np.array(["ACGT".index(base) for base in sequence], dtype=np.uint8)


# AGG at index 0; the only NGG in the string.
check("PAM at the very start is found at index 0",
      find_pam_starts(to_codes("AGGTTTT")).tolist(), [0])

# TGG starting at index 3 of "AAATGGAAA".
check("PAM in the middle is found at its own index",
      find_pam_starts(to_codes("AAATGGAAA")).tolist(), [3])

# "AAGGGAA" = A A G G G A A. Two NGG motifs overlap here:
#   index 1: A G G   (bases 1,2,3)
#   index 2: G G G   (bases 2,3,4)
check("GGG yields two overlapping PAM starts",
      find_pam_starts(to_codes("AAGGGAA")).tolist(), [1, 2])

print("\nTest 8: a guide is recovered from the sequence it came from")
# End-to-end: take a protospacer followed by a PAM, index it the way step 4
# does, and confirm the packed 20-mer equals the original protospacer.
protospacer = "GGAGAATGTCAGTCTGAGTC"
region = "TTTT" + protospacer + "AGG" + "TTTT"
codes = to_codes(region)
starts = find_pam_starts(codes)
weights = (4 ** np.arange(PROTOSPACER_LENGTH - 1, -1, -1)).astype(np.uint64)
recovered = []
for pam_start in starts:
    if pam_start >= PROTOSPACER_LENGTH:
        window = codes[pam_start - PROTOSPACER_LENGTH:pam_start]
        recovered.append(int(window.astype(np.uint64) @ weights))
check("the original protospacer is among the indexed sites",
      encode(protospacer) in recovered, True)

print()
if failures:
    print(f"{len(failures)} check(s) FAILED")
    sys.exit(1)
print("All checks passed.")
