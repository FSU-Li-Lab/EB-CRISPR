#!/usr/bin/env python3
"""
05_filter_qc.py

Apply QC thresholds and stage passing genomes for dereplication.

Output: 05_dereplication/input/  (dRep input)
        05_dereplication/high_quality_genomes.tsv
"""

import argparse
import shutil
from pathlib import Path

import pandas as pd


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qc", default="../04_qc/qc_summary.tsv")
    ap.add_argument("--fasta-dir", default="../03_genomes_filtered/fasta")
    ap.add_argument("--outdir", default="../05_dereplication/input")
    ap.add_argument("--min-completeness", type=float, default=95.0)
    ap.add_argument("--max-contamination", type=float, default=5.0)
    ap.add_argument("--symlink", action="store_true",
                    help="Symlink instead of copying (saves disk)")
    return ap.parse_args()


def main():
    args = parse_args()
    qc = pd.read_csv(args.qc, sep="\t")
    print(f"Genomes: {len(qc)}")

    good = qc[(qc.Completeness >= args.min_completeness) &
              (qc.Contamination <= args.max_contamination)]
    print(f"Passing QC: {len(good)}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    n = 0
    for genome in good.Genome:
        src = Path(args.fasta_dir) / f"{genome}.fna"
        dst = outdir / f"{genome}.fna"
        if not src.exists() or dst.exists():
            continue
        if args.symlink:
            dst.symlink_to(src.resolve())
        else:
            shutil.copy2(src, dst)
        n += 1

    good.to_csv(Path(args.outdir).parent / "high_quality_genomes.tsv",
                sep="\t", index=False)
    print(f"Staged: {n}")
    print(f"Output: {outdir}")
    print("\nNext: dRep dereplicate 05_dereplication/output \\")
    print("        -g 05_dereplication/input/*.fna -sa 0.99 -comp 95 -con 5 -p 80")


if __name__ == "__main__":
    main()
