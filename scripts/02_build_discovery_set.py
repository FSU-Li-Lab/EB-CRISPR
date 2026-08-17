#!/usr/bin/env python3
"""
02_build_discovery_set.py

Build the discovery set: all GTDB species representatives, plus a
genus-stratified sample of non-representatives.

Representatives are kept unconditionally so no species is lost. The stratified
sample on top gives within-species depth for the genera that matter clinically,
without letting Escherichia and Klebsiella - which dominate GenBank - crowd out
everything else.

NOTE ON GTDB GENUS SUFFIXES
---------------------------
GTDB splits polyphyletic genera with alphabetic suffixes: Klebsiella_A,
Escherichia_B, Citrobacter_A, Enterobacter_D and so on. Matching a genus by
exact string equality therefore silently skips every split lineage. This script
matches on the base genus name and includes all its suffixed variants, so
"Klebsiella" captures Klebsiella, Klebsiella_A, Klebsiella_B ...

Input:  01_metadata/Enterobacterales_high_quality_metadata.tsv
Output: 02_genomes_raw/discovery_set.tsv
        02_genomes_raw/discovery_accessions.txt

Usage:
    python3 02_build_discovery_set.py
    python3 02_build_discovery_set.py --seed 42
"""

import argparse
import re
from pathlib import Path

import pandas as pd

# Target counts per base genus, applied to NON-representatives.
# All GTDB species representatives are retained regardless.
TARGETS = {
    "Escherichia": 1200, "Klebsiella": 1200, "Salmonella": 1000,
    "Enterobacter": 700, "Citrobacter": 350, "Serratia": 300,
    "Yersinia": 300, "Proteus": 250, "Pantoea": 250,
    "Pectobacterium": 250, "Edwardsiella": 200, "Morganella": 150,
    "Providencia": 150, "Erwinia": 150, "Dickeya": 150, "Hafnia": 150,
}


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", default="../01_metadata/Enterobacterales_high_quality_metadata.tsv")
    ap.add_argument("--out-table", default="../02_genomes_raw/discovery_set.tsv")
    ap.add_argument("--out-accessions", default="../02_genomes_raw/discovery_accessions.txt")
    ap.add_argument("--seed", type=int, default=42)
    return ap.parse_args()


def base_genus(g):
    """Klebsiella_A -> Klebsiella"""
    return re.sub(r"_[A-Z]+$", "", g) if isinstance(g, str) else g


def main():
    args = parse_args()

    df = pd.read_csv(args.metadata, sep="\t", low_memory=False)
    df["genus"] = df["gtdb_taxonomy"].str.extract(r"g__([^;]*)")[0]
    df["base_genus"] = df["genus"].map(base_genus)
    print(f"High-quality genomes: {len(df)}")

    # All species representatives
    rep = df[df["gtdb_genome_representative"] == df["accession"]].copy()
    print(f"GTDB species representatives: {len(rep)}")

    # Genus-stratified sample of the rest
    rest = df[~df["accession"].isin(rep["accession"])]
    picked = []
    for genus, n in TARGETS.items():
        sub = rest[rest["base_genus"] == genus]
        if not len(sub):
            print(f"  {genus}: none available")
            continue
        take = min(n, len(sub))
        picked.append(sub.sample(n=take, random_state=args.seed))
        variants = sorted(sub["genus"].unique())
        print(f"  {genus}: {take} of {len(sub)}   ({', '.join(variants)})")

    discovery = pd.concat([rep] + picked).drop_duplicates("accession")
    print(f"\nDiscovery set: {len(discovery)} genomes")
    print(f"  species: {discovery['gtdb_taxonomy'].str.extract(r's__([^;]*)')[0].nunique()}")
    print(f"  genera:  {discovery['genus'].nunique()}")

    Path(args.out_table).parent.mkdir(parents=True, exist_ok=True)
    discovery.to_csv(args.out_table, sep="\t", index=False)

    acc = discovery["ncbi_genbank_assembly_accession"].dropna().drop_duplicates()
    acc.to_csv(args.out_accessions, index=False, header=False)

    print(f"\nSaved: {args.out_table}")
    print(f"       {args.out_accessions} ({len(acc)} accessions)")
    print("\nDownload with:")
    print(f"  datasets download genome accession --inputfile {args.out_accessions} \\")
    print("      --include genome --dehydrated --filename genomes.zip")
    print("  unzip -q genomes.zip -d raw && datasets rehydrate --directory raw")


if __name__ == "__main__":
    main()
