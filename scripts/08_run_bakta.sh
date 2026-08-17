#!/usr/bin/env bash
#
# 08_run_bakta.sh
#
# Annotate dRep representatives with Bakta.
#
# Usage:
#   DB=/path/to/bakta_db bash 06_run_bakta.sh
#   JOBS=4 THREADS=8 DB=/path/to/bakta_db bash 06_run_bakta.sh
#
set -euo pipefail

IN="${IN:-../05_dereplication/output/dereplicated_genomes}"
OUT="${OUT:-../06_annotation/bakta}"
DB="${DB:?set DB=/path/to/bakta_db}"
THREADS="${THREADS:-8}"
JOBS="${JOBS:-1}"

mkdir -p "$OUT"

run_one() {
    local genome="$1"
    local name
    name=$(basename "$genome" .fna)
    [[ -s "$OUT/$name/$name.tsv" ]] && return 0
    bakta --db "$DB" --output "$OUT/$name" --prefix "$name" \
          --threads "$THREADS" --force "$genome" \
          > "$OUT/$name.log" 2>&1 || echo "FAILED: $name" >> "$OUT/failed.txt"
}
export -f run_one
export OUT DB THREADS

echo "Genomes: $(ls -1 "$IN"/*.fna | wc -l)"

if command -v parallel >/dev/null 2>&1 && [[ "$JOBS" -gt 1 ]]; then
    ls -1 "$IN"/*.fna | parallel -j "$JOBS" --bar run_one {}
else
    for g in "$IN"/*.fna; do run_one "$g"; done
fi

echo "Done. Verify with: python3 09_check_bakta.py"
