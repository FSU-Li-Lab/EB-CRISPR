# 02_genomes_raw

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

Size: tens of GB. Gitignored.
