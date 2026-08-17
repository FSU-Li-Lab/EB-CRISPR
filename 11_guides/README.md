# 11_guides

Everything guide-related.

**Produced here:**
- `<locus>_candidates.tsv` / `.fasta` — step 12
- `all_candidates.fasta` — pooled with **locus-prefixed names**
- `candidates_nr.fasta` — after overlap collapse (step 16)
- `offtarget_tier1/`, `offtarget_full/` — screens (step 14)
- `site_counts*.tsv` — photon budget (step 18)
- `panel_final.fasta`, `order_sheet.tsv` — final output (step 20)

When pooling candidates, prefix names by locus. The scanner numbers guides
independently per file, so a plain `cat` silently collapses colliding names.
