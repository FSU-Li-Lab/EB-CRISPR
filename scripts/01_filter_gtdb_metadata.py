#!/usr/bin/env python3
"""
01_filter_gtdb_metadata.py

Select high-quality Enterobacterales genomes from the GTDB bacterial metadata
table.

GTDB is used for selection only - it distributes metadata, not assemblies.
The accessions produced downstream are fetched from NCBI (step 04).

Input:  01_metadata/bac120_metadata.tsv   (download from https://gtdb.ecogenomic.org/downloads)
Output: 01_metadata/Enterobacterales_high_quality_metadata.tsv

Usage:
    python3 01_filter_gtdb_metadata.py
    python3 01_filter_gtdb_metadata.py --min-completeness 95 --max-contamination 5
"""

import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", default="../01_metadata/bac120_metadata.tsv")
    ap.add_argument("--out", default="../01_metadata/Enterobacterales_high_quality_metadata.tsv")
    ap.add_argument("--taxon", default="o__Enterobacterales")
    ap.add_argument("--min-completeness", type=float, default=95.0)
    ap.add_argument("--max-contamination", type=float, default=5.0)
    return ap.parse_args()


def main():
    args = parse_args()

    print(f"Reading {args.metadata} ...")
    df = pd.read_csv(args.metadata, sep="\t", low_memory=False)
    print(f"  total GTDB genomes: {len(df)}")

    ent = df[df["gtdb_taxonomy"].str.contains(args.taxon, na=False)].copy()
    print(f"  {args.taxon}: {len(ent)}")

    hq = ent[
        (ent["checkm2_completeness"] >= args.min_completeness) &
        (ent["checkm2_contamination"] <= args.max_contamination)
    ].copy()
    print(f"  passing quality filter: {len(hq)}")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    hq.to_csv(args.out, sep="\t", index=False)
    print(f"\nSaved: {args.out}")


if __name__ == "__main__":
    main()
