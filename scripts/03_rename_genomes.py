#!/usr/bin/env python3
"""
03_rename_genomes.py

Flatten the NCBI datasets directory tree into one FASTA per accession.

Symlinks rather than copies, so the download can be re-verified against the
originals and the flattened tree costs no extra disk.

Input:  02_genomes_raw/raw/ncbi_dataset/data/<accession>/*_genomic.fna
Output: 03_genomes_filtered/fasta/<accession>.fna
        03_genomes_filtered/genome_mapping.tsv

Usage:
    python3 03_rename_genomes.py
"""

import argparse
import glob
import os
from pathlib import Path


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input-dir", default="../02_genomes_raw/raw/ncbi_dataset/data")
    ap.add_argument("--output-dir", default="../03_genomes_filtered/fasta")
    ap.add_argument("--mapping", default="../03_genomes_filtered/genome_mapping.tsv")
    return ap.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Exclude cds_from_genomic.fna and rna_from_genomic.fna - they match the
    # same glob but contain no rRNA / no intergenic sequence.
    files = [f for f in glob.glob(f"{args.input_dir}/**/*_genomic.fna", recursive=True)
             if not os.path.basename(f).startswith(("cds_from", "rna_from"))]
    print(f"Found genomes: {len(files)}")

    mapping = []
    for fasta in files:
        acc = os.path.basename(os.path.dirname(fasta))
        dst = outdir / f"{acc}.fna"
        if not dst.exists():
            os.symlink(os.path.abspath(fasta), dst)
        mapping.append((acc, os.path.abspath(fasta)))

    Path(args.mapping).parent.mkdir(parents=True, exist_ok=True)
    with open(args.mapping, "w") as fh:
        fh.write("Genome_ID\tOriginal_path\n")
        for acc, path in mapping:
            fh.write(f"{acc}\t{path}\n")

    print(f"Linked: {len(mapping)}")
    print(f"Output: {outdir}")
    print(f"Mapping: {args.mapping}")


if __name__ == "__main__":
    main()
