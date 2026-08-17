#!/usr/bin/env bash
#
# setup_dirs.sh
#
# Create the working directory tree the pipeline expects.
#
# Scripts run from scripts/ and use ../ relative paths, so every data directory
# is a sibling of scripts/. Each gets a README explaining what belongs in it and
# which step produces it. All are gitignored - only the structure is tracked.
#
# Usage:
#   bash setup_dirs.sh
#
set -euo pipefail

ROOT="${1:-.}"
cd "$ROOT"

mk() {
    mkdir -p "$1"
    printf '%s\n' "$2" > "$1/README.md"
    touch "$1/.gitkeep"
}

mk 01_metadata '# 01_metadata

**You provide:** `bac120_metadata.tsv` — download from
https://gtdb.ecogenomic.org/downloads (the bacterial metadata table, ~2 GB
uncompressed). Nothing else is needed to start.

**Produced here:**
- `Enterobacterales_high_quality_metadata.tsv` — step 01

GTDB is used for *selection only*. It distributes metadata, not assemblies;
genomes come from NCBI in step 02.'

mk 02_genomes_raw '# 02_genomes_raw

**Produced here:**
- `discovery_set.tsv`, `discovery_accessions.txt` — step 02
- `raw/` — the NCBI datasets download tree

Download after step 02:

```bash
datasets download genome accession --inputfile discovery_accessions.txt \
    --include genome --dehydrated --filename genomes.zip
unzip -q genomes.zip -d raw
datasets rehydrate --directory raw --max-workers 10
```

Always dehydrated. A single large archive is one long HTTP/2 stream that dies
partway and leaves an unusable zip; rehydration is resumable.

Size: tens of GB. Gitignored.'

mk 03_genomes_filtered '# 03_genomes_filtered

**Produced here:**
- `fasta/` — one symlink per accession, `<accession>.fna` (step 03)
- `genome_mapping.tsv` — accession to original path

Symlinks, not copies, so this costs no extra disk. Step 03 excludes
`cds_from_genomic.fna` and `rna_from_genomic.fna`, which match the same glob
but contain no rRNA.'

mk 04_qc '# 04_qc

**You produce** by running the tools:
- `checkm2/quality_report.tsv`
- `quast/report.tsv`

```bash
checkm2 predict --input ../03_genomes_filtered/fasta --output-directory checkm2 \
    -x fna --threads 40
quast.py ../03_genomes_filtered/fasta/*.fna -o quast --threads 40
```

**Produced here:** `qc_summary.tsv` — step 04

Note: CheckM2 completeness comes from single-copy *protein* markers and is
blind to rRNA. A genome can score 100% complete with every rRNA operon
collapsed. For rRNA-targeted design, assembly level matters more.'

mk 05_dereplication '# 05_dereplication

**Produced here:**
- `input/` — QC-passing genomes staged for dRep (step 05)
- `high_quality_genomes.tsv` — step 05
- `output/dereplicated_genomes/` — representatives (step 06)
- `representatives_composition.tsv` — audit (step 07)

`output/dereplicated_genomes/` is the input to everything downstream.

dRep at 0.99 ANI is the long pole in the pipeline — budget many hours.'

mk 06_annotation '# 06_annotation

**You provide:** a Bakta database (~30 GB), path passed as `DB=`.

```bash
bakta_db download --output /path/to/bakta_db --type full
```

**Produced here:**
- `bakta/<genome>/` — one directory per genome (step 08)
- `summary/bakta_status.tsv`, `summary/bakta_rerun.txt` — step 09

Used downstream for `.tsv` (gene names) and `.ffn` (nucleotide CDS) in step 11.
rRNA guides do not depend on the annotation.'

mk 09_background '# 09_background

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

Size: >100 GB at default settings. Gitignored.'

mk 10_rrna '# 10_rrna

**Produced here (step 10):**
- `per_genome/` — barrnap output per genome
- `all_rrna.fa` — pooled
- `16S_full.fa`, `23S_full.fa` — length-filtered
- `16S_one_per_genome.fa`, `23S_one_per_genome.fa` — **input to step 12**

One copy per genome, since the operons within a genome are near-identical and
keeping all ~7 would weight that genome sevenfold in the alignment.

Expect to lose genomes here — short-read assemblies truncate rRNA operons and
the length filter drops them. Check retention *per genus* before proceeding.'

mk 10_core_genes '# 10_core_genes

**Produced here (step 11):** one FASTA per gene — `gyrB.fna`, `rpoB.fna`,
`secY.fna`, `atpD.fna`, ...

Nucleotide sequence, not protein: a guide is 20 nt of DNA, and amino-acid
conservation says nothing about whether any DNA window is conserved. In
practice most housekeeping genes yield zero order-wide candidates for exactly
this reason.

Keep only genes recovered from >=95% of genomes.'

mk 11_guides '# 11_guides

Everything guide-related.

**Produced here:**
- `<locus>_candidates.tsv` / `.fasta` — step 12
- `all_candidates.fasta` — pooled with **locus-prefixed names**
- `candidates_nr.fasta` — after overlap collapse (step 16)
- `offtarget_tier1/`, `offtarget_full/` — screens (step 14)
- `site_counts*.tsv` — photon budget (step 18)
- `panel_final.fasta`, `order_sheet.tsv` — final output (step 20)

When pooling candidates, prefix names by locus. The scanner numbers guides
independently per file, so a plain `cat` silently collapses colliding names.'

mkdir -p logs
printf '%s\n' '# logs

Run logs and status files. Gitignored.' > logs/README.md
touch logs/.gitkeep

echo "Created:"
for d in 01_metadata 02_genomes_raw 03_genomes_filtered 04_qc 05_dereplication \
         06_annotation 09_background 10_rrna 10_core_genes 11_guides logs; do
    echo "  $d/"
done
echo
echo "Next: put bac120_metadata.tsv in 01_metadata/, then"
echo "  cd scripts && python3 01_filter_gtdb_metadata.py"
