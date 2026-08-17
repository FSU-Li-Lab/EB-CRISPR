# 10_rrna

**Produced here (step 10):**
- `per_genome/` — barrnap output per genome
- `all_rrna.fa` — pooled
- `16S_full.fa`, `23S_full.fa` — length-filtered
- `16S_one_per_genome.fa`, `23S_one_per_genome.fa` — **input to step 12**

One copy per genome, since the operons within a genome are near-identical and
keeping all ~7 would weight that genome sevenfold in the alignment.

Expect to lose genomes here — short-read assemblies truncate rRNA operons and
the length filter drops them. Check retention *per genus* before proceeding.
