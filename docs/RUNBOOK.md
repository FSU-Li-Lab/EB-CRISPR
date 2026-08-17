# Runbook

Full command sequence with checkpoints. All scripts run **from `scripts/`** and
use `../` relative paths, so set up a working tree alongside the repo:

```
project/
├── scripts/          <- this repo's scripts/
├── 01_metadata/
├── 02_genomes_raw/
├── 03_genomes_filtered/
├── 04_qc/
├── 05_dereplication/
├── 06_annotation/
├── 09_background/
├── 10_rrna/
├── 10_core_genes/
└── 11_guides/
```

```bash
conda env create -f environment.yml && conda activate grna-design
bash setup_dirs.sh          # create the data directory tree
conda create -n casoff -c conda-forge -c bioconda cas-offinder
cd scripts
```

---

## 1. Dataset construction

Download `bac120_metadata.tsv` from https://gtdb.ecogenomic.org/downloads into
`01_metadata/`.

```bash
python3 01_filter_gtdb_metadata.py
python3 02_build_discovery_set.py
```

**Checkpoint:** the printed genus list should show suffixed variants
(`Klebsiella, Klebsiella_A, ...`). If only base names appear, suffix handling is
broken and split lineages are being skipped.

Download **dehydrated** — a single large archive is one long HTTP/2 stream that
frequently dies partway and leaves an unusable zip. Rehydration fetches files
individually and is resumable.

```bash
cd ../02_genomes_raw
datasets download genome accession --inputfile discovery_accessions.txt \
    --include genome --dehydrated --filename genomes.zip
unzip -q genomes.zip -d raw
datasets rehydrate --directory raw --max-workers 10

find raw -name "*_genomic.fna" | wc -l    # compare to accession count
cd ../scripts
python3 03_rename_genomes.py
```

Some accessions may be suppressed at NCBI and will never rehydrate. Record the
count; it belongs in your methods.

Runtime: hours, network-bound.

---

## 2. QC and dereplication

Run CheckM2 and QUAST over `03_genomes_filtered/fasta/`, then:

```bash
python3 04_build_qc_summary.py
python3 05_filter_qc.py --symlink
THREADS=80 bash 06_dereplicate.sh
python3 07_check_composition.py
```

**Checkpoint:** `07_check_composition.py` reports species and genera retained —
the dereplication audit a reviewer will ask for. If a genus vanished,
investigate before proceeding.

Runtime: dRep is the long pole, many hours to a day.

---

## 3. Annotation

```bash
DB=/path/to/bakta_db JOBS=4 THREADS=8 bash 08_run_bakta.sh
python3 09_check_bakta.py
```

**Checkpoint:** `Incomplete` and `Never started` must both be 0. The script also
flags genomes whose `.faa` is far below median — usually a reduced-genome
symbiont or a truncated assembly.

---

## 4. Target loci

Independent; run in either order.

```bash
JOBS=20 THREADS=2 bash 10_extract_rrna.sh
python3 11_extract_core_genes.py
```

**Checkpoint:** rRNA retention should exceed ~90% of genomes, core genes ≥95%.
Anything lower means the alignment under-represents part of the clade — check
*per genus* before continuing.

Runtime: 1–2 h for rRNA, ~15 min for core genes.

---

## 5. Guide discovery

Stratify by **genus**. With GTDB collapsing an entire order into one family,
family stratification enumerates from a handful of genomes and misses most of
the candidate space.

```bash
mkdir -p ../11_guides

for gene in 16S 23S; do
    python3 12_scan_conserved_grna.py \
        --fasta ../10_rrna/${gene}_one_per_genome.fa \
        --out ../11_guides/${gene}_candidates.tsv \
        --min-coverage 0.95 --refs-rank genus --refs-per-family 3
done

for gene in gyrB rpoB infB atpD recA secY; do
    [[ -s ../10_core_genes/${gene}.fna ]] || continue
    python3 12_scan_conserved_grna.py \
        --fasta ../10_core_genes/${gene}.fna \
        --out ../11_guides/${gene}_candidates.tsv \
        --min-coverage 0.80 --refs-rank genus
done
```

Pool with **locus-prefixed names** — the scanner numbers guides independently
per file, so a plain `cat` silently collapses colliding names:

```bash
cd ../11_guides
for f in *_candidates.fasta; do
    locus=$(basename "$f" _candidates.fasta)
    awk -v L="$locus" '/^>/{sub(/^>/,""); print ">"L"_"$1; next}{print}' "$f"
done > all_candidates.fasta
grep -c '^>' all_candidates.fasta
cd ../scripts

python3 16_check_guide_overlap.py --guides ../11_guides/all_candidates.fasta \
    --out ../11_guides/candidates_nr.fasta
```

Collapsing overlaps before screening typically removes 60–70% of candidates and
saves proportional compute.

**Checkpoint:** zero core-gene candidates is a real biological result, not a
failure — conserved proteins often contain no conserved 20 nt DNA window.

---

## 6. Background construction

```bash
python3 13_build_background_set.py --tiers 1 2 --cap-per-order 150
```

Genomes are capped per GTDB order rather than sampled at random, because risk
tracks phylogenetic proximity and gut abundance — neither of which correlates
with how many genomes an order happens to have deposited. Orders adjacent to
the target and the dominant gut orders are exempt from the cap.

**Checkpoint:** sister orders (Pasteurellales, Vibrionales) and Bacteroidales
must appear near the top of the order table. If environmental orders dominate,
the cap did not apply.

Download dehydrated as in step 1, then flatten with **genus in the contig
header** — off-target hits cannot be abundance-weighted afterwards otherwise,
and retrofitting means re-screening:

```bash
cd ../09_background
mkdir -p fasta

python3 - <<'PY'
import pandas as pd, re
m = pd.read_csv("background_metadata.tsv", sep="\t", low_memory=False)
with open("acc_genus.tsv", "w") as out:
    for _, r in m.iterrows():
        acc = r.get("ncbi_genbank_assembly_accession")
        if isinstance(acc, str):
            g = re.search(r"g__([^;]*)", str(r.get("gtdb_taxonomy", "")))
            out.write(f"{acc}\t{g.group(1) if g else 'unknown'}\n")
PY

while IFS=$'\t' read -r acc genus; do
    f=$(find raw/ncbi_dataset/data/"$acc" -name "*_genomic.fna" 2>/dev/null | head -1)
    [[ -n "$f" ]] || continue
    awk -v g="$genus" -v a="$acc" \
        '/^>/{n++; print ">"g"|"a"|"n; next}{print}' "$f" > fasta/"$acc".fna
done < acc_genus.tsv
cd ../scripts
```

Optionally add the host genome, flattened the same way into
`09_background/human/` with `>Homo|GCF_...|n` headers.

---

## 7. Off-target screen

Stage it. Screen sister orders first — small, fast, highest risk — then take
survivors to the full background.

```bash
conda activate casoff

mkdir -p ../09_background/tier1_fasta
while read acc; do
    [[ -f ../09_background/fasta/$acc.fna ]] && \
        ln -sf ../fasta/$acc.fna ../09_background/tier1_fasta/
done < ../09_background/background_tier1_accessions.txt

GUIDES=../11_guides/candidates_nr.fasta \
BG=../09_background/tier1_fasta \
OUT=../11_guides/offtarget_tier1 \
JOBS=20 MM_STRICT=2 VALIDATE=0 bash 14_offtarget_screen.sh

conda activate grna-design
python3 15_summarise_offtargets.py \
    --offtarget-dir ../11_guides/offtarget_tier1 \
    --guides ../11_guides/candidates_nr.fasta --gate 10
```

Then the full background on survivors:

```bash
conda activate casoff
GUIDES=../11_guides/offtarget_tier1/surviving_candidates.fasta \
BG=../09_background/fasta \
OUT=../11_guides/offtarget_full \
JOBS=20 MM_STRICT=2 VALIDATE=0 bash 14_offtarget_screen.sh
```

Use `JOBS=1` for a large host genome — cas-offinder loads the whole sequence
into memory, and several copies at once will exhaust it.

**Checkpoint:** verify the run actually happened. Compare
`find OUT/strict -name '*.txt' | wc -l` against the genome count, and check a
log for "Comparing patterns". An empty input directory produces "0 hits", which
is indistinguishable from a clean result.

Runtime: ~30 min per 32k genomes for a handful of guides at 20 jobs.

---

## 8. Panel assembly

```bash
conda activate grna-design

python3 18_count_sites.py \
    --guides ../11_guides/offtarget_full/surviving_candidates.fasta \
    --out ../11_guides/site_counts.tsv --closed-only --gate 10 --threads 40
```

`--closed-only` restricts to `Complete Genome` assemblies. Essential for
multicopy loci: scaffolded records carry N-gaps exactly where rRNA repeats sit,
and CheckM2 cannot detect this because it scores protein markers.

Rank by per-genus minimum coverage, not by mean — a guide at 99% overall can be
entirely absent from one genus:

```bash
python3 - <<'PY'
import pandas as pd
d = pd.read_csv('../11_guides/site_counts.tsv', sep='\t')
d['genus'] = d.species.str.split().str[0]
g = [c for c in d.columns if c.startswith(('16S_','23S_','secY_','gRNA_'))]
cov = d.groupby('genus')[g].apply(lambda x: (x > 0).mean())
print([c for c in g if (cov[c] == 1.0).all()])
PY
```

Build the panel from guides at full coverage across every in-scope genus, then
validate and generate sequences:

```bash
python3 16_check_guide_overlap.py --guides ../11_guides/panel_final.fasta
python3 17_check_spacing.py --guides ../11_guides/panel_final.fasta \
    --locus 23S ../10_rrna/23S_one_per_genome.fa \
    --locus secY ../10_core_genes/secY.fna
python3 18_count_sites.py --guides ../11_guides/panel_final.fasta \
    --out ../11_guides/site_counts_FINAL.tsv --closed-only --gate 10 --threads 40
python3 20_make_order_sheet.py --guides ../11_guides/panel_final.fasta \
    --out ../11_guides/order_sheet.tsv
```

`19_select_panel.py` automates selection under an abundance-weighted purity
budget when the candidate set is too large to inspect by hand.

**Before ordering:** confirm each protospacer is followed by NGG in a reference
sequence, and independently verify the sgRNA scaffold sequence against its
primary source.
