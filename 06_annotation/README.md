# 06_annotation

**You provide:** a Bakta database (~30 GB), path passed as `DB=`.

```bash
bakta_db download --output /path/to/bakta_db --type full
```

**Produced here:**
- `bakta/<genome>/` — one directory per genome (step 08)
- `summary/bakta_status.tsv`, `summary/bakta_rerun.txt` — step 09

Used downstream for `.tsv` (gene names) and `.ffn` (nucleotide CDS) in step 11.
rRNA guides do not depend on the annotation.
