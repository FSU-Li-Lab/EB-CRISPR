#!/usr/bin/env python3
"""
11_extract_core_genes.py

Pull single-copy housekeeping loci out of the Bakta annotation you already
have. These become the specificity half of the two-class panel (proposal
section 2.4): one site per genome each, but discriminating cleanly at family
and order level, where rRNA guides discriminate poorly.

Runs alongside 21_extract_rrna.sh, not after it. rRNA supplies photons, core
genes supply specificity, and 25_select_panel.py trades them off.

Reads {genome}.tsv for the locus-tag -> gene-name mapping and {genome}.ffn for
the nucleotide CDS. Nucleotide, not protein: a guide is a 20mer of DNA, and
the amino-acid conservation of rpoB tells you nothing about whether any 20 nt
window inside it is conserved.

Output: one FASTA per gene, headers as {genome}|{locus_tag}, ready for
12_scan_conserved_grna.py.

Usage:
    python 11_extract_core_genes.py
    python 11_extract_core_genes.py --genes gyrB rpoB infB atpD recA
"""

import argparse
import re
from pathlib import Path

# Housekeeping loci named in the proposal, plus the Enterobacteriaceae MLST
# scheme genes and a few universally single-copy backups.
DEFAULT_GENES = [
    "gyrB", "rpoB", "infB", "atpD",            # proposal section 2.4
    "adk", "fumC", "icd", "mdh", "purA", "recA",   # E. coli MLST scheme
    "secY", "groL", "dnaK", "rpoD", "tuf", "fusA",
]

RE_SUFFIX = re.compile(r"_\d+$")


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bakta-dir", default="../06_annotation/bakta")
    ap.add_argument("--outdir", default="../10_core_genes")
    ap.add_argument("--genes", nargs="+", default=DEFAULT_GENES)
    ap.add_argument("--min-len", type=int, default=300,
                    help="Drop obvious fragments")
    return ap.parse_args()


def read_bakta_tsv_genes(path, wanted):
    """locus_tag -> canonical gene name, for the genes we care about."""
    header, out = None, {}
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if header is None:
                if line.startswith("#") and "Sequence" in line and "Type" in line:
                    header = [c.lstrip("#").strip() for c in line.split("\t")]
                continue
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < len(header):
                f += [""] * (len(header) - len(f))
            row = dict(zip(header, f))
            gene = (row.get("Gene") or "").strip()
            if not gene:
                continue
            # Bakta suffixes paralogues: rpoB_1, rpoB_2
            base = RE_SUFFIX.sub("", gene)
            key = wanted.get(base.lower())
            if key:
                lt = (row.get("Locus Tag") or "").strip()
                if lt:
                    out[lt] = key
    return out


def read_fasta(path):
    name, seq = None, []
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq).upper()
                name, seq = line[1:], []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq).upper()


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    wanted = {g.lower(): g for g in args.genes}
    handles = {g: open(outdir / f"{g}.fna", "w") for g in args.genes}
    counts = {g: 0 for g in args.genes}
    multi = {g: 0 for g in args.genes}

    dirs = sorted(p for p in Path(args.bakta_dir).iterdir() if p.is_dir())
    print(f"Bakta annotations: {len(dirs)}")
    print(f"Genes: {', '.join(args.genes)}\n")

    n_ok = 0
    for n, d in enumerate(dirs, 1):
        name = d.name
        tsv, ffn = d / f"{name}.tsv", d / f"{name}.ffn"
        if not (tsv.exists() and ffn.exists()):
            continue

        lt_to_gene = read_bakta_tsv_genes(tsv, wanted)
        if not lt_to_gene:
            continue
        n_ok += 1

        seen = set()
        for header, seq in read_fasta(ffn):
            lt = header.split()[0]
            gene = lt_to_gene.get(lt)
            if not gene or len(seq) < args.min_len:
                continue
            # One copy per genome. A second hit means a paralogue or a
            # duplicated region; taking the first keeps the alignment
            # one-row-per-genome, which is what the scanner expects.
            if gene in seen:
                multi[gene] += 1
                continue
            seen.add(gene)
            handles[gene].write(f">{name}|{lt}\n{seq}\n")
            counts[gene] += 1

        if n % 500 == 0:
            print(f"  {n}/{len(dirs)}")

    for h in handles.values():
        h.close()

    print(f"\nGenomes with usable annotation: {n_ok}\n")
    print(f"{'gene':<8} {'genomes':>8} {'%':>7} {'extra copies':>13}")
    for g in args.genes:
        pct = 100 * counts[g] / n_ok if n_ok else 0
        print(f"{g:<8} {counts[g]:>8} {pct:>6.1f}% {multi[g]:>13}")

    print(f"\nSaved to {outdir}")
    print("\nA gene below ~95% is not a usable single-copy target - either the")
    print("annotation is inconsistent across families or the gene genuinely")
    print("varies. Feed the high-recovery ones to 12_scan_conserved_grna.py.")


if __name__ == "__main__":
    main()
