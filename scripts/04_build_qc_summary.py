#!/usr/bin/env python3
"""
04_build_qc_summary.py

Merge CheckM2 and QUAST reports into one QC table.

IMPORTANT: CheckM2 completeness is computed from single-copy PROTEIN markers,
so it is structurally blind to rRNA operons. A genome can score 100% complete
while every rRNA repeat has been collapsed or replaced by an assembly gap. For
rRNA-targeted guide design, contig count and assembly level matter more than
CheckM2 completeness - see the note in 18_count_sites.py.

Input:  04_qc/checkm2/quality_report.tsv, 04_qc/quast/report.tsv
Output: 04_qc/qc_summary.tsv
"""

import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--checkm2", default="../04_qc/checkm2/quality_report.tsv")
    ap.add_argument("--quast", default="../04_qc/quast/report.tsv")
    ap.add_argument("--out", default="../04_qc/qc_summary.tsv")
    return ap.parse_args()


def main():
    args = parse_args()

    checkm = pd.read_csv(args.checkm2, sep="\t")
    checkm = checkm.rename(columns={"Name": "Genome"})
    keep = [c for c in ["Genome", "Completeness", "Contamination", "Genome_Size",
                        "GC_Content", "Total_Contigs", "Contig_N50"] if c in checkm.columns]
    qc = checkm[keep].copy()

    if Path(args.quast).exists():
        quast = pd.read_csv(args.quast, sep="\t").T
        quast.columns = quast.iloc[0]
        quast = quast.iloc[1:]
        quast.index.name = "Genome"
        quast = quast.reset_index()

        rename = {}
        for c in quast.columns:
            if "# contigs" in str(c):
                rename[c] = "QUAST_Contigs"
            elif "Total length" in str(c):
                rename[c] = "QUAST_Length"
            elif "N50" == str(c):
                rename[c] = "QUAST_N50"
        quast = quast.rename(columns=rename)
        cols = ["Genome"] + [c for c in ["QUAST_Contigs", "QUAST_Length", "QUAST_N50"]
                             if c in quast.columns]
        qc = qc.merge(quast[cols], on="Genome", how="left")
    else:
        print(f"QUAST report not found at {args.quast}; using CheckM2 only")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    qc.to_csv(args.out, sep="\t", index=False)
    print(f"Genomes: {len(qc)}")
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
