#!/usr/bin/env bash
#
# 14_offtarget_screen.sh
#
# Replacement for the single-shot version when cas-offinder segfaults.
#
# cas-offinder loads the whole background into memory before searching. Pointed
# at a directory of several thousand genomes with a few hundred query patterns,
# it will either exhaust RAM or die on the first malformed file, and in both
# cases the failure looks identical: a segfault with no useful message. Running
# one genome per invocation costs a little startup time and removes both
# problems. It also isolates a bad file to a single line in failed.txt instead
# of killing the run, and it parallelises cleanly.
#
# Usage:
#   GUIDES=../11_guides/all_candidates.fasta \
#   BG=../09_background/fasta \
#   OUT=../11_guides/offtarget_tier1 \
#   JOBS=20 bash 23b_offtarget_batched.sh
#
set -euo pipefail

ROOT="${ROOT:-..}"
GUIDES="${GUIDES:-$ROOT/11_guides/all_candidates.fasta}"
BG="${BG:-$ROOT/09_background/fasta}"
OUT="${OUT:-$ROOT/11_guides/offtarget}"
DEVICE="${DEVICE:-C}"
JOBS="${JOBS:-20}"
MM_STRICT="${MM_STRICT:-6}"
VALIDATE="${VALIDATE:-1}"

mkdir -p "$OUT/strict" "$OUT/seed" "$OUT/logs"

command -v cas-offinder >/dev/null 2>&1 || {
    echo "ERROR: cas-offinder not on PATH" >&2; exit 1; }
[[ -s "$GUIDES" ]] || { echo "ERROR: $GUIDES missing or empty" >&2; exit 1; }

# ----------------------------------------------------------------------
# Pre-flight validation
#
# Incomplete rehydration is the usual culprit: `datasets rehydrate` leaves
# zero-byte placeholders for files it did not fetch, and cas-offinder
# segfaults on them rather than reporting an error.
# ----------------------------------------------------------------------
if [[ "$VALIDATE" == "1" ]]; then
    echo "Validating background files ..."
    : > "$OUT/logs/bad_files.txt"
    nbad=0
    while IFS= read -r f; do
        real=$(readlink -f "$f" 2>/dev/null || echo "$f")
        if [[ ! -s "$real" ]]; then
            echo -e "$f\tempty_or_broken" >> "$OUT/logs/bad_files.txt"; nbad=$((nbad+1)); continue
        fi
        if ! head -c 1 "$real" | grep -q '>'; then
            echo -e "$f\tnot_fasta" >> "$OUT/logs/bad_files.txt"; nbad=$((nbad+1)); continue
        fi
        if ! grep -q '[ACGTNacgtn]' "$real"; then
            echo -e "$f\tno_sequence" >> "$OUT/logs/bad_files.txt"; nbad=$((nbad+1)); continue
        fi
    done < <(find "$BG" -maxdepth 1 -name "*.fna")
    echo "  unusable files: $nbad  (see $OUT/logs/bad_files.txt)"
    if (( nbad > 0 )); then
        echo "  These are excluded. If the count is large, rehydration is"
        echo "  incomplete - rerun 'datasets rehydrate --directory raw'."
    fi
fi

# ----------------------------------------------------------------------
# Query files, built once and reused by every invocation
# ----------------------------------------------------------------------
grep -v '^>' "$GUIDES" | sed 's/[[:space:]]*$//' | grep -v '^$' > "$OUT/guides.txt"
n_guides=$(wc -l < "$OUT/guides.txt")
echo "Guides: $n_guides"

: > "$OUT/queries_strict.txt"
: > "$OUT/queries_seed.txt"
while read -r s; do
    echo -e "${s}NNN\t${MM_STRICT}" >> "$OUT/queries_strict.txt"
    seed="${s: -12}"
    echo -e "NNNNNNNN${seed}NNN\t0" >> "$OUT/queries_seed.txt"
done < "$OUT/guides.txt"

run_one() {
    local genome="$1"
    local name mode qfile outfile inp
    name=$(basename "$genome" .fna)

    for mode in strict; do
        outfile="$OUT/$mode/$name.txt"
        [[ -f "$outfile" ]] && continue
        qfile="$OUT/queries_${mode}.txt"
        inp=$(mktemp)
        {
            readlink -f "$genome"
            echo "NNNNNNNNNNNNNNNNNNNNNRG"
            cat "$qfile"
        } > "$inp"

        if cas-offinder "$inp" "$DEVICE" "$outfile" \
                > "$OUT/logs/$name.$mode.log" 2>&1; then
            :
        else
            echo -e "$name\t$mode" >> "$OUT/logs/failed.txt"
            rm -f "$outfile"
        fi
        rm -f "$inp"
    done
}
export -f run_one
export OUT DEVICE

mapfile -t files < <(find "$BG" -maxdepth 1 -name "*.fna" \
    | grep -v -F -f <(cut -f1 "$OUT/logs/bad_files.txt" 2>/dev/null || true) 2>/dev/null \
    || find "$BG" -maxdepth 1 -name "*.fna")

echo "Background genomes: ${#files[@]}"
echo "Jobs: $JOBS"
echo

printf '%s\n' "${files[@]}" | parallel -j "$JOBS" --bar run_one {}

# ----------------------------------------------------------------------
# Pool
# ----------------------------------------------------------------------
echo
for mode in strict; do
    find "$OUT/$mode" -name "*.txt" -exec cat {} + > "$OUT/offtargets_${mode}.txt" 2>/dev/null || : > "$OUT/offtargets_${mode}.txt"
    echo "$mode hits: $(wc -l < "$OUT/offtargets_${mode}.txt")"
done

if [[ -f "$OUT/logs/failed.txt" ]]; then
    echo
    echo "Failed genomes: $(wc -l < "$OUT/logs/failed.txt")  (see $OUT/logs/failed.txt)"
    echo "Rerun this script to retry only those - completed outputs are skipped."
fi

echo
echo "Saved to $OUT"
