#!/usr/bin/env python3
"""
07_check_composition.py

Report the taxonomic composition of the dereplicated representative set.

This is the dereplication audit a reviewer will ask for: how many species and
genera survived, and whether any clade was lost. Run it before committing CPU
to annotation.

Output: 05_dereplication/representatives_composition.tsv
"""

import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", default="../05_dereplication/output/dereplicated_genomes")
    ap.add_argument("--metadata", default="../02_genomes_raw/discovery_set.tsv")
    ap.add_argument("--out", default="../05_dereplication/representatives_composition.tsv")
    return ap.parse_args()


def main():
    args = parse_args()

    reps = sorted(p.stem for p in Path(args.genomes).glob("*.fna"))
    m = pd.read_csv(args.metadata, sep="\t", low_memory=False)
    m = m.rename(columns={"ncbi_genbank_assembly_accession": "Genome"}).drop_duplicates("Genome")

    d = pd.DataFrame({"Genome": reps}).merge(m, on="Genome", how="left")
    d["genus"] = d.gtdb_taxonomy.str.extract(r"g__([^;]*)")
    d["species"] = d.gtdb_taxonomy.str.extract(r"s__([^;]*)")

    print(f"Representatives:  {len(d)}")
    print(f"Distinct species: {d.species.nunique()}")
    print(f"Distinct genera:  {d.genus.nunique()}")

    if "ncbi_assembly_level" in d.columns:
        print(f"\nAssembly level:\n{d.ncbi_assembly_level.value_counts().to_string()}")
    if "checkm2_completeness" in d.columns:
        print(f"\nMean completeness:  {d.checkm2_completeness.mean():.1f}%")
        print(f"Mean contamination: {d.checkm2_contamination.mean():.2f}%")

    print(f"\nGenus distribution:\n{d.genus.value_counts().head(25).to_string()}")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(args.out, sep="\t", index=False)
    print(f"\nSaved: {args.out}")


if __name__ == "__main__":
    main()
