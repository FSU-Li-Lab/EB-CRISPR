#!/usr/bin/env python3
"""
09_check_bakta.py

Verify that Bakta actually finished for every genome.

The Bakta loop in 06_run_bakta.sh runs with --force and redirects stderr into a
per-genome log, so a genome that crashed (OOM, truncated contig, DB hiccup)
leaves behind a directory that *looks* fine. This script checks the output
files instead of the directory names, and writes a re-run list for anything
incomplete.

Usage:
    python 09_check_bakta.py
    python 09_check_bakta.py --bakta-dir ../06_annotation/bakta
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

# Files Bakta writes for every genome. If any is missing or empty, the run
# did not complete cleanly.
REQUIRED_SUFFIXES = [".tsv", ".faa", ".ffn", ".gff3", ".txt"]


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bakta-dir", default="../06_annotation/bakta",
                    help="Directory holding one sub-directory per genome")
    ap.add_argument("--genome-dir", default="../05_dereplication/output/dereplicated_genomes",
                    help="dRep representative FASTAs, used to find genomes Bakta never started")
    ap.add_argument("--outdir", default="../06_annotation/summary")
    return ap.parse_args()


def main():
    args = parse_args()

    bakta_dir = Path(args.bakta_dir)
    genome_dir = Path(args.genome_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not bakta_dir.is_dir():
        sys.exit(f"ERROR: {bakta_dir} does not exist")

    # ------------------------------------------------------------------
    # Expected genomes
    # ------------------------------------------------------------------
    expected = sorted(p.stem for p in genome_dir.glob("*.fna")) if genome_dir.is_dir() else []
    if expected:
        print(f"Genomes in {genome_dir}: {len(expected)}")
    else:
        print(f"WARNING: no FASTAs found in {genome_dir} - "
              f"cannot detect genomes Bakta never started")

    sample_dirs = sorted(p for p in bakta_dir.iterdir() if p.is_dir())
    print(f"Bakta output directories:  {len(sample_dirs)}")
    print()

    # ------------------------------------------------------------------
    # Check each output directory
    # ------------------------------------------------------------------
    rows = []
    for d in sample_dirs:
        name = d.name
        missing = []
        sizes = {}

        for suffix in REQUIRED_SUFFIXES:
            f = d / f"{name}{suffix}"
            if not f.exists():
                missing.append(suffix)
            elif f.stat().st_size == 0:
                missing.append(suffix + "(empty)")
            else:
                sizes[suffix] = f.stat().st_size

        rows.append({
            "Genome": name,
            "Status": "OK" if not missing else "INCOMPLETE",
            "Missing": ";".join(missing),
            "faa_bytes": sizes.get(".faa", 0),
            "tsv_bytes": sizes.get(".tsv", 0),
        })

    df = pd.DataFrame(rows)

    # ------------------------------------------------------------------
    # Genomes that never produced a directory at all
    # ------------------------------------------------------------------
    never_started = []
    if expected:
        have = set(df["Genome"]) if len(df) else set()
        never_started = [g for g in expected if g not in have]

    # ------------------------------------------------------------------
    # Report
    # ------------------------------------------------------------------
    n_ok = int((df["Status"] == "OK").sum()) if len(df) else 0
    n_bad = int((df["Status"] == "INCOMPLETE").sum()) if len(df) else 0

    print(f"Complete:      {n_ok}")
    print(f"Incomplete:    {n_bad}")
    print(f"Never started: {len(never_started)}")
    print()

    if n_bad:
        print("Incomplete genomes (first 20):")
        print(df[df["Status"] == "INCOMPLETE"][["Genome", "Missing"]].head(20).to_string(index=False))
        print()

    if never_started:
        print("Never started (first 20):")
        for g in never_started[:20]:
            print("  " + g)
        print()

    # Flag suspiciously small .faa files - usually a truncated or
    # near-empty assembly that slipped through QC.
    if n_ok:
        ok = df[df["Status"] == "OK"]
        cutoff = ok["faa_bytes"].median() * 0.25
        tiny = ok[ok["faa_bytes"] < cutoff]
        if len(tiny):
            print(f"WARNING: {len(tiny)} genomes have a .faa < 25% of the median size. "
                  f"Inspect these before annotation analysis:")
            print(tiny[["Genome", "faa_bytes"]].head(20).to_string(index=False))
            print()

    # ------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------
    status_file = outdir / "bakta_status.tsv"
    df.to_csv(status_file, sep="\t", index=False)

    rerun = list(df[df["Status"] == "INCOMPLETE"]["Genome"]) + never_started
    rerun_file = outdir / "bakta_rerun.txt"
    with open(rerun_file, "w") as fh:
        for g in rerun:
            fh.write(g + "\n")

    print("Saved:")
    print(" ", status_file)
    print(" ", rerun_file, f"({len(rerun)} genomes)")

    if rerun:
        print()
        print("Re-run with:")
        print(f"  while read g; do bash 06_run_bakta_one.sh "
              f"{genome_dir}/$g.fna; done < {rerun_file}")


if __name__ == "__main__":
    main()
