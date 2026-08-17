# 09_background

The negative set. Specificity is defined entirely by what is here.

**Produced here:**
- `background_metadata.tsv`, `background_accessions.txt` — step 13
- `background_tier1_accessions.txt` — sister orders, screen these first
- `raw/` — NCBI download tree
- `fasta/` — flattened, **with genus in the contig header**
- `tier1_fasta/` — symlinks for the staged first-pass screen
- `human/GRCh38.fna` — optional host genome

Contig headers must be `>genus|accession|n`. Off-target hits cannot be
abundance-weighted without the genus, and retrofitting means re-screening.
The flattening loop is in `docs/RUNBOOK.md` step 6.

Size: >100 GB at default settings. Gitignored.
