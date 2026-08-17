# 04_qc

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
collapsed. For rRNA-targeted design, assembly level matters more.
