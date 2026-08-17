# 03_genomes_filtered

**Produced here:**
- `fasta/` — one symlink per accession, `<accession>.fna` (step 03)
- `genome_mapping.tsv` — accession to original path

Symlinks, not copies, so this costs no extra disk. Step 03 excludes
`cds_from_genomic.fna` and `rna_from_genomic.fna`, which match the same glob
but contain no rRNA.
