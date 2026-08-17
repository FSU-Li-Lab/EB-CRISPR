#!/usr/bin/env bash
#
# 10_extract_rrna.sh
#
# Extract 16S and 23S rRNA from every genome, with a length filter.
#
# Why rRNA, and why a length filter matters here:
#
# rRNA operons are the only loci in a bacterial genome that are conserved at
# the NUCLEOTIDE level across an entire order while still carrying
# order-discriminating regions. They are also multi-copy (7 operons in E. coli
# K-12), which for an amplification-free dCas9 biosensor is worth roughly an
# order of magnitude in sensitivity. That combination is why the guide panel
# will almost certainly land here.
#
# The catch: short-read assemblies collapse and truncate rRNA operons. rRNA
# repeats are longer than the read insert, so assemblers break contigs inside
# them. A large fraction of GenBank draft genomes carry a partial 16S, or one
# split across two contigs, or none at all. If you feed those into an alignment
# you will conclude that a perfectly conserved region is variable, and discard
# your best guide. Hence --min-16s / --min-23s: partial genes are dropped, not
# repaired.
#
# Expect to lose genomes here. That is correct behaviour. Check the retention
# rate per family before proceeding - if one family drops out entirely, your
# order-level coverage claim is not supported by the data.
#
# Install:  mamba install -c bioconda barrnap seqkit
# Usage:    bash 21_extract_rrna.sh
#
set -euo pipefail

ROOT="${ROOT:-..}"
IN="${IN:-$ROOT/05_dereplication/output/dereplicated_genomes}"
OUT="${OUT:-$ROOT/10_rrna}"
JOBS="${JOBS:-20}"
THREADS="${THREADS:-2}"

# Full-length 16S is ~1540 bp, 23S ~2900 bp. Allow some slack for
# barrnap's boundary calls but reject obvious fragments.
MIN_16S="${MIN_16S:-1400}"
MIN_23S="${MIN_23S:-2600}"

mkdir -p "$OUT/per_genome" "$OUT/logs"

for tool in barrnap seqkit; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: $tool not on PATH" >&2; exit 1; }
done

run_one() {
    local genome="$1"
    local name
    name=$(basename "$genome" .fna)

    [[ -f "$OUT/per_genome/$name.fa" ]] && return 0

    barrnap --kingdom bac --threads "$THREADS" --quiet \
            --outseq "$OUT/per_genome/$name.tmp.fa" \
            "$genome" > "$OUT/per_genome/$name.gff" 2>"$OUT/logs/$name.log" || {
        echo "FAILED: $name" >> "$OUT/logs/failed.txt"; return 0; }

    # barrnap headers look like: >16S_rRNA::contig_1:12345-13884(+)
    # Prefix the genome id so sequences stay traceable after pooling.
    if [[ -s "$OUT/per_genome/$name.tmp.fa" ]]; then
        sed "s/^>/>${name}|/" "$OUT/per_genome/$name.tmp.fa" \
            > "$OUT/per_genome/$name.fa"
    else
        : > "$OUT/per_genome/$name.fa"
    fi
    rm -f "$OUT/per_genome/$name.tmp.fa"
}
export -f run_one
export OUT THREADS

n_total=$(ls -1 "$IN"/*.fna 2>/dev/null | wc -l)
echo "Genomes: $n_total"
echo

if command -v parallel >/dev/null 2>&1; then
    ls -1 "$IN"/*.fna | parallel -j "$JOBS" --bar run_one {}
else
    ls -1 "$IN"/*.fna | xargs -P "$JOBS" -I{} bash -c 'run_one "$@"' _ {}
fi

echo
echo "Pooling and filtering ..."

cat "$OUT"/per_genome/*.fa > "$OUT/all_rrna.fa"

# 16S
seqkit grep -nrp "16S_rRNA" "$OUT/all_rrna.fa" \
    | seqkit seq -m "$MIN_16S" > "$OUT/16S_full.fa"

# 23S
seqkit grep -nrp "23S_rRNA" "$OUT/all_rrna.fa" \
    | seqkit seq -m "$MIN_23S" > "$OUT/23S_full.fa"

# One representative copy per genome. The operons within a genome are near
# identical, and keeping all 7 would weight that genome 7x in the alignment.
for gene in 16S 23S; do
    seqkit fx2tab "$OUT/${gene}_full.fa" \
        | awk -F'\t' '{split($1,a,"|"); if (!(a[1] in seen)) {seen[a[1]]=1; print}}' \
        | seqkit tab2fx > "$OUT/${gene}_one_per_genome.fa"
done

echo
echo "Retention:"
for gene in 16S 23S; do
    total=$(grep -c '^>' "$OUT/${gene}_full.fa" || true)
    uniq=$(grep -c '^>' "$OUT/${gene}_one_per_genome.fa" || true)
    echo "  $gene: $total copies, $uniq genomes of $n_total ($((100 * uniq / n_total))%)"
done

echo
echo "Saved to $OUT"
echo
echo "IMPORTANT: check retention per family before continuing --"
echo "  grep '^>' $OUT/16S_one_per_genome.fa | cut -d'|' -f1 | sed 's/^>//' \\"
echo "    > $OUT/genomes_with_16S.txt"
echo "then join against your GTDB taxonomy. A family at low retention means"
echo "your alignment underrepresents it, not that it lacks 16S."
echo
echo "Next: align, then run 12_scan_conserved_grna.py"
echo "  mafft --auto --thread $JOBS $OUT/16S_one_per_genome.fa > $OUT/16S_aln.fa"
