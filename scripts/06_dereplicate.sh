#!/usr/bin/env bash
#
# 06_dereplicate.sh
#
# Dereplicate at 0.99 ANI with dRep.
#
# Why 0.99 and not lower: it collapses near-identical clinical isolates so that
# over-sequenced lineages (E. coli ST131, K. pneumoniae ST258) do not dominate
# conservation estimates, while retaining genuine within-species diversity.
#
# Caveat worth knowing: dereplication removes strains you still intend to
# detect. Final guide performance is therefore measured on the FULL QC-passing
# set, not on the dereplicated representatives (see 18_count_sites.py).
#
# Usage:
#   bash 07_dereplicate.sh
#   THREADS=80 bash 07_dereplicate.sh
#
set -euo pipefail

IN="${IN:-../05_dereplication/input}"
OUT="${OUT:-../05_dereplication/output}"
ANI="${ANI:-0.99}"
THREADS="${THREADS:-40}"

command -v dRep >/dev/null 2>&1 || { echo "ERROR: dRep not on PATH" >&2; exit 1; }

echo "Input genomes: $(ls -1 "$IN"/*.fna | wc -l)"
echo "ANI threshold: $ANI"

dRep dereplicate "$OUT" \
    -g "$IN"/*.fna \
    -sa "$ANI" \
    -comp 95 \
    -con 5 \
    -p "$THREADS"

echo
echo "Representatives: $(ls -1 "$OUT/dereplicated_genomes"/*.fna | wc -l)"
