"""Step 3 — Apply standard sequence filters to the candidate guides.

None of these filters predict how well a guide will cut. They remove guides
that are likely to fail for known mechanical reasons, which is a different
thing. Real on-target scoring (Rule Set 2/3, Doench et al. 2016) is a trained
model; see the README for why this project does not reimplement it.

Filters applied
---------------
GC content 40-80%
    Very low GC binds the target weakly; very high GC increases off-target
    binding and can be hard to synthesise. 40-80% is the usual window.

No TTTT run
    Guides are normally transcribed from a U6 promoter, and a run of four
    or more T's is a termination signal for RNA polymerase III. A TTTT in
    the protospacer can cut the guide RNA short.

No homopolymer of 5 or more
    Long single-base runs are associated with synthesis errors and poorer
    activity.

Overlaps a coding exon
    For a knockout you want to disrupt protein-coding sequence, so guides
    in introns or UTRs are dropped here.

Output:
  results/filtered_guides.tsv   guides that passed every filter
  results/filter_summary.txt    how many guides each filter removed
  results/guide_positions.png   where the surviving guides sit in the gene
"""

import re
import tomllib
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # no interactive window; we only save a file
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO

ROOT = Path(__file__).resolve().parent.parent
config = tomllib.loads((ROOT / "config.toml").read_text())
gene = config["gene"]

guides = pd.read_csv(ROOT / "results" / "all_sites.tsv", sep="\t")

# Coding exons, for the figure. Transcripts overlap heavily, so merge the
# union of their exons into non-overlapping blocks.
record = SeqIO.read(ROOT / "data" / f"{gene}.gb", "genbank")
exon_parts = sorted(
    (int(part.start), int(part.end))
    for feature in record.features if feature.type == "CDS"
    for part in feature.location.parts
)
coding_exons = []
for start, end in exon_parts:
    if coding_exons and start <= coding_exons[-1][1]:
        coding_exons[-1][1] = max(coding_exons[-1][1], end)
    else:
        coding_exons.append([start, end])


def gc_percent(protospacer: str) -> float:
    return 100 * (protospacer.count("G") + protospacer.count("C")) / len(protospacer)


def longest_homopolymer(protospacer: str) -> int:
    return max(len(run) for run in re.findall(r"(A+|C+|G+|T+)", protospacer))


guides["gc_percent"] = guides["protospacer"].map(gc_percent)
guides["longest_run"] = guides["protospacer"].map(longest_homopolymer)
guides["has_tttt"] = guides["protospacer"].str.contains("TTTT")

checks = {
    "in a coding exon": guides["in_cds"],
    "GC between 40% and 80%": guides["gc_percent"].between(40, 80),
    "no TTTT terminator": ~guides["has_tttt"],
    "no run of 5+ identical bases": guides["longest_run"] < 5,
}

# Report each filter's effect on its own, then apply them together.
summary_lines = [f"Starting candidates: {len(guides):,}", ""]
for description, passed in checks.items():
    summary_lines.append(
        f"  {description:<32} keeps {passed.sum():>6,}  "
        f"(removes {(~passed).sum():>6,})"
    )

keep = pd.Series(True, index=guides.index)
for passed in checks.values():
    keep &= passed

filtered = guides[keep].copy()
# A guide that cuts near the start of the coding sequence is more likely to
# produce a frameshift that knocks the protein out, so sort by cut position.
filtered.sort_values("cut_site", inplace=True)

summary_lines += ["", f"Passed every filter:  {len(filtered):,}"]
summary = "\n".join(summary_lines)

results = ROOT / "results"
filtered.to_csv(results / "filtered_guides.tsv", sep="\t", index=False)
(results / "filter_summary.txt").write_text(summary + "\n", encoding="utf-8")

# --- Figure: where the surviving guides sit --------------------------------
figure, (axis_top, axis_bottom) = plt.subplots(
    2, 1, figsize=(11, 5), height_ratios=[1, 1.4], constrained_layout=True
)

sequence_length = int(guides["end"].max())

# Bottom lane: the gene body with its coding exons marked.
axis_top.hlines(-0.5, 0, sequence_length, color="#d8d8d8", linewidth=9)
for exon_start, exon_end in coding_exons:
    axis_top.hlines(-0.5, exon_start, exon_end, color="#3d7ea6", linewidth=9)

# Upper lane: one tick per surviving guide, split by strand.
for _, guide in filtered.iterrows():
    if guide["strand"] == "+":
        axis_top.vlines(guide["cut_site"], 0.15, 0.85, color="#2a9d5c", linewidth=0.8)
    else:
        axis_top.vlines(guide["cut_site"], -0.15, 0.1, color="#d97b3a", linewidth=0.8)

axis_top.set_title(
    f"{gene}: {len(filtered):,} guides passing all filters "
    "(green = + strand, orange = - strand; blue = coding exons)"
)
axis_top.set_xlim(0, sequence_length)
axis_top.set_ylim(-0.9, 1.1)
axis_top.set_yticks([])
axis_top.set_xlabel("position in the downloaded gene region (bp)")

axis_bottom.hist(guides["gc_percent"], bins=40, color="#cccccc",
                 label="all candidates")
axis_bottom.hist(filtered["gc_percent"], bins=40, color="#2a9d5c",
                 label="passed filters")
axis_bottom.axvline(40, color="#b03030", linestyle="--", linewidth=1)
axis_bottom.axvline(80, color="#b03030", linestyle="--", linewidth=1)
axis_bottom.set_xlabel("GC content of the 20-nt protospacer (%)")
axis_bottom.set_ylabel("number of guides")
axis_bottom.legend()

figure.savefig(results / "guide_positions.png", dpi=150)

print(summary)
print("\nSaved: results/filtered_guides.tsv, results/filter_summary.txt, "
      "results/guide_positions.png")
