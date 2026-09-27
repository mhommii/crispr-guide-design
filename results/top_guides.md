# Top 10 candidate guides for TP53

Ranked by predicted specificity and cut position. **On-target cutting efficiency is not predicted here** - see the note in `steps/05_rank_and_report.py`.

| Rank | Protospacer (5'->3') | PAM | Strand | GC% | Cut position in CDS | Off-targets 1-2 mm | 3-4 mm | Intact seed |
|---|---|---|---|---|---|---|---|---|
| 1 | `TCTCGAAGCGCTCACGCCCA` | CGG | + | 65 | 16% | 0 | 0 | 0 |
| 2 | `CCACGGATCTGCAGCAACAG` | AGG | + | 60 | 16% | 0 | 3 | 0 |
| 3 | `CCTGAGTAGTGGTAATCTAC` | TGG | - | 45 | 60% | 0 | 1 | 0 |
| 4 | `AGATTACCACTACTCAGGAT` | AGG | + | 40 | 61% | 0 | 1 | 0 |
| 5 | `CAGCCACCTGAAGTCCAAAA` | AGG | - | 50 | 1% | 0 | 4 | 0 |
| 6 | `GCCTGTGTTATCTCCTAGGT` | TGG | - | 50 | 67% | 0 | 1 | 0 |
| 7 | `GCCAACCTAGGAGATAACAC` | AGG | + | 50 | 67% | 0 | 1 | 0 |
| 8 | `TGCACGGTCAGTTGCCCTGA` | GGG | - | 60 | 91% | 0 | 0 | 0 |
| 9 | `GGGCAGCTACGGTTTCCGTC` | TGG | - | 65 | 92% | 0 | 0 | 0 |
| 10 | `CCCCGGACGATATTGAACAA` | TGG | - | 50 | 95% | 0 | 0 | 0 |

Searched against chromosome 17 only (223 guides passed the step 3 filters). Counts exclude each guide's own target site.
