"""Step 5 — Rank the guides and write the final report.

What this ranking is, and what it is not
----------------------------------------
This ranks guides on *specificity* - how likely they are to cut somewhere
they should not - plus position within the gene. It deliberately does not
predict *on-target activity*, i.e. how efficiently a guide cuts its
intended target. Activity prediction needs a model trained on large screens
(Rule Set 2, Doench et al. 2016; CRISPick uses this). Reimplementing such a
model from memory would produce numbers that look authoritative and are
not, so it is left out. For real work, take the shortlist from here and put
the sequences through CRISPick.

The published specificity scores (the MIT score and the CFD score) each use
a matrix of experimentally measured mismatch penalties. Those matrices are
not reproduced here either. Instead the ranking below uses only what this
project actually measured, and says exactly how each part is weighted.

Ranking, in order of importance
-------------------------------
1. Off-targets with an intact seed. The 12 bases next to the PAM must pair
   for Cas9 to cut, so an off-target whose seed matches perfectly is the
   most likely to be cut in reality. Weighted most heavily.
2. Off-targets at 1-2 mismatches, which are also plausible cut sites.
3. Off-targets at 3-4 mismatches, which are much less likely; counted, but
   weighted lightly.
4. Position in the coding sequence. A cut early in the protein-coding
   region is more likely to knock the protein out through a frameshift,
   so earlier cuts rank higher as a tie-breaker.

Output:
  results/ranked_guides.tsv   all guides, best first
  results/top_guides.md       the top 10 as a readable table
  results/offtarget_profile.png  specificity of the top guides
"""

import tomllib
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO

ROOT = Path(__file__).resolve().parent.parent
config = tomllib.loads((ROOT / "config.toml").read_text())
gene = config["gene"]
results = ROOT / "results"

guides = pd.read_csv(results / "filtered_guides.tsv", sep="\t")
summary = pd.read_csv(results / "offtarget_summary.tsv", sep="\t")
guides = guides.merge(summary, on="protospacer", how="left")

# Where does the coding sequence start? Used to score cut position.
record = SeqIO.read(ROOT / "data" / f"{gene}.gb", "genbank")
coding_starts = [int(part.start)
                 for feature in record.features if feature.type == "CDS"
                 for part in feature.location.parts]
coding_start, coding_end = min(coding_starts), max(
    int(part.end) for feature in record.features if feature.type == "CDS"
    for part in feature.location.parts
)
coding_span = coding_end - coding_start

# Fraction of the way through the coding region, 0 = start, 1 = end.
guides["position_in_cds"] = (
    (guides["cut_site"] - coding_start) / coding_span
).clip(0, 1)

# --- Penalty score: lower is better ----------------------------------------
# Weights are choices, not measurements. They are written here explicitly so
# anyone reading can disagree with them and change them.
WEIGHT_INTACT_SEED = 40   # off-target with a perfectly matching seed
WEIGHT_CLOSE = 10         # off-target at 1-2 mismatches
WEIGHT_DISTANT = 1        # off-target at 3-4 mismatches
WEIGHT_POSITION = 5       # how late in the protein the cut falls

guides["penalty"] = (
    WEIGHT_INTACT_SEED * guides["offtargets_intact_seed"]
    + WEIGHT_CLOSE * (guides["mm1"] + guides["mm2"])
    + WEIGHT_DISTANT * (guides["mm3"] + guides["mm4"])
    + WEIGHT_POSITION * guides["position_in_cds"]
)

guides.sort_values(["penalty", "cut_site"], inplace=True)
guides.reset_index(drop=True, inplace=True)
guides.insert(0, "rank", guides.index + 1)
guides.to_csv(results / "ranked_guides.tsv", sep="\t", index=False)

# --- Readable top-10 table -------------------------------------------------
top = guides.head(10)
lines = [
    f"# Top 10 candidate guides for {gene}",
    "",
    "Ranked by predicted specificity and cut position. **On-target cutting "
    "efficiency is not predicted here** - see the note in "
    "`steps/05_rank_and_report.py`.",
    "",
    "| Rank | Protospacer (5'->3') | PAM | Strand | GC% | Cut position in CDS "
    "| Off-targets 1-2 mm | 3-4 mm | Intact seed |",
    "|---|---|---|---|---|---|---|---|---|",
]
for _, guide in top.iterrows():
    lines.append(
        f"| {guide['rank']} | `{guide['protospacer']}` | {guide['pam']} | "
        f"{guide['strand']} | {guide['gc_percent']:.0f} | "
        f"{guide['position_in_cds']:.0%} | "
        f"{int(guide['mm1'] + guide['mm2'])} | "
        f"{int(guide['mm3'] + guide['mm4'])} | "
        f"{int(guide['offtargets_intact_seed'])} |"
    )
lines += [
    "",
    f"Searched against chromosome 17 only ({len(guides)} guides passed the "
    "step 3 filters). Counts exclude each guide's own target site.",
]
(results / "top_guides.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

# --- Figure: specificity profile of the top guides -------------------------
figure, axis = plt.subplots(figsize=(10, 5), constrained_layout=True)
display = guides.head(20).iloc[::-1]
labels = [f"{row['rank']:>2}. {row['protospacer']}" for _, row in display.iterrows()]
positions = range(len(display))

axis.barh(positions, display["mm1"] + display["mm2"], color="#b03030",
          label="off-targets, 1-2 mismatches")
axis.barh(positions, display["mm3"] + display["mm4"],
          left=display["mm1"] + display["mm2"], color="#e0b080",
          label="off-targets, 3-4 mismatches")
axis.set_yticks(list(positions))
axis.set_yticklabels(labels, fontfamily="monospace", fontsize=8)
axis.set_xlabel("off-target sites found on chromosome 17")
axis.set_title(f"{gene}: specificity of the 20 best-ranked guides")
# Put the legend outside the plotting area so it cannot cover a bar.
axis.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)
axis.margins(x=0.02)
figure.savefig(results / "offtarget_profile.png", dpi=150)

print(f"Ranked {len(guides)} guides.")
print(f"Coding sequence spans {coding_start:,}-{coding_end:,} "
      f"in the downloaded region.")
print(f"\nBest guide: {guides.iloc[0]['protospacer']} "
      f"({guides.iloc[0]['pam']}, {guides.iloc[0]['strand']} strand)")
print(f"  off-targets 1-2 mismatches: "
      f"{int(guides.iloc[0]['mm1'] + guides.iloc[0]['mm2'])}")
print(f"  off-targets 3-4 mismatches: "
      f"{int(guides.iloc[0]['mm3'] + guides.iloc[0]['mm4'])}")
print(f"  cut at {guides.iloc[0]['position_in_cds']:.0%} through the CDS")
print("\nSaved: results/ranked_guides.tsv, results/top_guides.md, "
      "results/offtarget_profile.png")
