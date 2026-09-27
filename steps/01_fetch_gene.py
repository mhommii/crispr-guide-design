"""Step 1 — Download a gene's genomic sequence and exon map from NCBI.

Why genomic DNA and not mRNA?
Cas9 cuts DNA in the genome. A guide copied from mRNA could span an
exon–exon junction that doesn't exist in the chromosome, so it would never
bind. We therefore download the gene's region of the chromosome (GRCh38),
together with the annotation that says where the coding exons are.

Output:
  data/<GENE>.gb     GenBank record: sequence + feature annotations
  data/<GENE>.fasta  the same sequence in plain FASTA format
"""

import sys
import tomllib
from pathlib import Path

from Bio import Entrez, SeqIO

ROOT = Path(__file__).resolve().parent.parent
config = tomllib.loads((ROOT / "config.toml").read_text())

if "@" not in config["email"]:
    sys.exit("Put your email in config.toml first (NCBI requires it).")

Entrez.email = config["email"]
gene, organism = config["gene"], config["organism"]

# 1. Find the NCBI Gene record, e.g. TP53 -> Gene ID 7157
query = f'{gene}[Gene Name] AND "{organism}"[Organism] AND alive[property]'
with Entrez.esearch(db="gene", term=query) as handle:
    ids = Entrez.read(handle)["IdList"]
if not ids:
    sys.exit(f"No NCBI Gene record found for {gene} in {organism}.")
gene_id = ids[0]

# 2. Ask where the gene sits on the reference chromosome
with Entrez.esummary(db="gene", id=gene_id) as handle:
    summary = Entrez.read(handle)["DocumentSummarySet"]["DocumentSummary"][0]
loc = summary["GenomicInfo"][0]
chrom_acc = loc["ChrAccVer"]  # e.g. NC_000017.11 = chromosome 17, GRCh38
strand = "-" if int(loc["ChrStart"]) > int(loc["ChrStop"]) else "+"
start, stop = sorted((int(loc["ChrStart"]), int(loc["ChrStop"])))
# NCBI gives 0-based coordinates; efetch expects 1-based, inclusive.
seq_start, seq_stop = start + 1, stop + 1

# 3. Download that region with its annotations
with Entrez.efetch(db="nuccore", id=chrom_acc, rettype="gbwithparts",
                   retmode="text", seq_start=seq_start, seq_stop=seq_stop) as handle:
    record = SeqIO.read(handle, "genbank")

data_dir = ROOT / "data"
data_dir.mkdir(exist_ok=True)
SeqIO.write(record, data_dir / f"{gene}.gb", "genbank")
SeqIO.write(record, data_dir / f"{gene}.fasta", "fasta")

cds_features = [f for f in record.features if f.type == "CDS"]
print(f"Gene:        {gene} (NCBI Gene ID {gene_id}) — {summary['Description']}")
print(f"Location:    {chrom_acc}:{seq_start:,}-{seq_stop:,} (gene on {strand} strand)")
print(f"Length:      {len(record.seq):,} bp")
print(f"CDS records: {len(cds_features)} (one per annotated transcript/isoform)")
print(f"Saved:       data/{gene}.gb, data/{gene}.fasta")
