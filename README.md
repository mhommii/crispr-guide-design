# CRISPR Guide Design

Designing and checking SpCas9 guide RNAs for a human cancer-associated gene (default: **TP53**), one documented step at a time.

> **How this was made:** a guided learning project built and run with Claude Code (an AI assistant), which wrote the code and the explanations and executed each step. It is a learning exercise, not independent research. The *My notes* sections are left for me to fill in as I work through the code.

**Status:** all 5 steps complete for TP53.

## What it does

```text
gene symbol → genomic sequence + exons → NGG PAM sites → candidate sgRNAs
            → sequence filters → off-target search → ranked guide table
```

| Step | Script | What it does |
|---|---|---|
| 1 | `steps/01_fetch_gene.py` | Downloads the gene's GRCh38 region and exon annotation from NCBI |
| 2 | `steps/02_find_pam_sites.py` | Finds every 20-nt protospacer next to an NGG PAM, on both strands |
| 3 | `steps/03_filter_guides.py` | Filters on GC content, TTTT terminator, homopolymers, coding-exon overlap |
| 4 | `steps/04_off_target_search.py` | Searches chromosome 17 for sites the guide could also cut |
| 5 | `steps/05_rank_and_report.py` | Ranks by specificity and cut position, writes the report |

## Results for TP53

| | |
|---|---|
| Region downloaded | `chr17:7,668,421-7,687,490` (19,070 bp, minus-strand gene) |
| NGG sites found | 2,860 (1,391 on +, 1,469 on −) |
| Passed all filters | 223 (280 were in coding exons before the other filters) |
| Sites searched on chr17 | 10,495,968 NGG sites, both strands |
| Guides with no off-target at 1–2 mismatches | 205 of 223 |
| Guides with no off-target at all (≤4 mismatches) | 7 |

**Top-ranked guide:** `TCTCGAAGCGCTCACGCCCA` (PAM `CGG`, + strand, 65% GC, cuts 16% into the coding sequence) — no off-target sites at up to 4 mismatches anywhere on chromosome 17.

Full table: [`results/top_guides.md`](results/top_guides.md) · all 223 ranked: [`results/ranked_guides.tsv`](results/ranked_guides.tsv)

![Guide positions across TP53](results/guide_positions.png)

Guides cluster exactly on the coding exons (blue), which is the coding-exon filter working as intended.

![Off-target profile of top guides](results/offtarget_profile.png)

## Checks

The two error-prone parts have known-answer tests, where the expected result can be worked out by hand:

```powershell
.venv\Scripts\python.exe tests\test_pam_logic.py          # PAM finding, both strands
.venv\Scripts\python.exe tests\test_mismatch_counting.py  # packed mismatch counting
```

`test_mismatch_counting.py` also compares the fast bit-packed comparison against plain base-by-base counting on 20,000 random sequence pairs.

Step 4 additionally self-checks: every guide came from TP53 on chromosome 17, so every guide must find its own target site there. It stops with an error if any guide does not. That check caught a real off-by-one bug in the PAM indexing during development — before the fix, 154 of 223 guides silently failed to find themselves.

## Running it

```powershell
py -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
# put your own email in config.toml (NCBI requires one for downloads)

.venv\Scripts\python.exe steps\01_fetch_gene.py
.venv\Scripts\python.exe steps\02_find_pam_sites.py
.venv\Scripts\python.exe steps\03_filter_guides.py
# step 4 needs the reference chromosome (~25 MB download, 81 MB unpacked):
#   https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr17.fa.gz -> data/
.venv\Scripts\python.exe steps\04_off_target_search.py   # ~90 s
.venv\Scripts\python.exe steps\05_rank_and_report.py
```

To target a different gene, change `gene` in `config.toml`. Note that step 4's reference chromosome must be the one that gene sits on.

## My notes

*(to be written by me)*

### Step 1 — fetching the gene

### Step 2 — PAM sites

### Step 3 — filters

### Step 4 — off-targets

### Step 5 — ranking

## Limitations

These matter, and the numbers above should not be read without them.

- **Learning project.** The guides are **not experimentally validated** and are not recommendations for any therapeutic use.
- **On-target activity is not predicted.** How *well* a guide cuts is a separate question needing a model trained on large screens (Rule Set 2; used by CRISPick). This project ranks specificity and position only. The published MIT and CFD specificity scores are also not reproduced, because both depend on experimentally measured mismatch penalty matrices. The ranking weights used instead are stated explicitly in `steps/05_rank_and_report.py` — they are reasonable choices, not measurements.
- **One chromosome, not the genome.** Off-target search covers chr17 (~2.7% of the genome). A guide that looks clean here has not been checked against the other 97%.
- **No bulges.** Only mismatches are searched, not insertions or deletions between guide and DNA. Cas-OFFinder can search those.
- **Cas-OFFinder was not run here.** It is the standard tool and step 4 writes an input file for it (`results/cas_offinder_input.txt`), but it needs an OpenCL runtime that this machine does not have. Step 4 implements the same mismatch search directly instead — exhaustive over chr17, so not an approximation, but also not a validated reimplementation of that tool.
- **Annotation-dependent.** "Coding exon" means the union of all annotated CDS features in the downloaded RefSeq record. A guide in an exon used by only some transcripts is treated the same as one in a constitutive exon.

## References

- Cock, P. J. A. et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/btp163
- Bae, S., Park, J. & Kim, J.-S. (2014). Cas-OFFinder: a fast and versatile algorithm that searches for potential off-target sites of Cas9 RNA-guided endonucleases. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/btu048 · [snugel/cas-offinder](https://github.com/snugel/cas-offinder)
- Doench, J. G. et al. (2016). Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. *Nature Biotechnology*. https://doi.org/10.1038/nbt.3437
- Sanson, K. R. et al. (2018). Optimized libraries for CRISPR-Cas9 genetic screens with multiple modalities. *Nature Communications*. https://doi.org/10.1038/s41467-018-07901-8
- CRISPick (Broad Institute Genetic Perturbation Platform), web tool: https://portals.broadinstitute.org/gppx/crispick/public
- Sequence data: NCBI Gene and RefSeq (https://www.ncbi.nlm.nih.gov/gene/); reference chromosome from the UCSC Genome Browser hg38 release (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/)
