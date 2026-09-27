"""Known-answer test for the PAM-finding logic in step 2.

Run:  .venv\\Scripts\\python.exe tests\\test_pam_logic.py

The point is to check the rules by hand on a sequence small enough to read,
before trusting the same code on 19 kb of TP53.
"""

import re
import sys

from Bio.Seq import Seq

PROTOSPACER_LENGTH = 20


def find_sites(sequence: str):
    """Same logic as steps/02_find_pam_sites.py, kept standalone for testing."""
    sequence = sequence.upper()
    sites = []
    for match in re.finditer(r"(?=([ACGT]GG))", sequence):
        start = match.start() - PROTOSPACER_LENGTH
        if start < 0:
            continue
        sites.append(("+", sequence[start:match.start()], match.group(1), start))
    for match in re.finditer(r"(?=(CC[ACGT]))", sequence):
        end = match.start() + 3 + PROTOSPACER_LENGTH
        if end > len(sequence):
            continue
        forward = sequence[match.start() + 3:end]
        sites.append(("-", str(Seq(forward).reverse_complement()),
                      str(Seq(match.group(1)).reverse_complement()),
                      match.start() + 3))
    return sites


failures = []


def check(description, got, expected):
    if got == expected:
        print(f"  PASS  {description}")
    else:
        print(f"  FAIL  {description}\n        expected {expected!r}\n        got      {got!r}")
        failures.append(description)


print("Test 1: one PAM on the forward strand")
# 20 A's, then AGG. The only NGG here is the AGG at position 20, so the
# protospacer must be the 20 A's.
seq1 = "A" * 20 + "AGG"
sites1 = find_sites(seq1)
forward1 = [s for s in sites1 if s[0] == "+"]
check("finds exactly one forward site", len(forward1), 1)
check("protospacer is the 20 preceding bases", forward1[0][1], "A" * 20)
check("PAM is AGG", forward1[0][2], "AGG")
check("protospacer starts at position 0", forward1[0][3], 0)

print("\nTest 2: a PAM too close to the start is skipped")
# Only 5 bases before the AGG: not enough room for a 20-nt protospacer.
check("no site when fewer than 20 bases precede the PAM",
      [s for s in find_sites("ACGTA" + "AGG") if s[0] == "+"], [])

print("\nTest 3: overlapping PAMs are both reported")
# ...GGG contains AGG (at index 20) and GGG (at index 21) -> two sites.
seq3 = "C" * 20 + "AGGG"
check("GGG yields two overlapping forward sites",
      len([s for s in find_sites(seq3) if s[0] == "+"]), 2)

print("\nTest 4: reverse strand")
# CCT at the start reads as AGG on the reverse strand. The 20 bases after it
# (all A on the forward strand) become 20 T's once reverse complemented.
seq4 = "CCT" + "A" * 20
reverse4 = [s for s in find_sites(seq4) if s[0] == "-"]
check("finds exactly one reverse site", len(reverse4), 1)
check("protospacer is reverse complemented", reverse4[0][1], "T" * 20)
check("PAM reads AGG on the reverse strand", reverse4[0][2], "AGG")

print("\nTest 5: the guide really pairs with the target")
# Independent check: take the reverse-strand hit and confirm that its
# reverse complement is present in the original sequence.
check("reverse complement of the guide is found in the input",
      str(Seq(reverse4[0][1]).reverse_complement()) in seq4, True)

print("\nTest 6: a sequence with no NGG produces nothing")
check("no sites in a G-free sequence", find_sites("ACTACTACTACTACTACTACTACT"), [])

print()
if failures:
    print(f"{len(failures)} check(s) FAILED")
    sys.exit(1)
print("All checks passed.")
