#!/usr/bin/env python3
"""
24_count_sites.py

Count dCas9 binding SITES per genome for each candidate guide.

For a detection assay, presence/absence of a target is what matters. For
Hy-SCALE it is not: the sorter reads photons, and photons scale with the
number of labeled sites per bead. So the quantity that determines whether a
cell crosses the FACS gate is

    S(genome) = sum over guides of (sites for that guide in that genome)

which is the [S] bracket in the photon budget of the proposal, §2.5.

This distinction changes what a "good" guide is. A guide present once in every
genome in the order is worse for enrichment than a guide present seven times in
90% of them, even though the first looks better by prevalence. Presence is
necessary; copy number is what buys signal.

The second consequence is the one that governs M1.2. Abundance is preserved
only if every target species sits comfortably above the gate. A species whose
S is near threshold is partially recovered, and partial recovery of one species
and full recovery of another is exactly the compositional distortion the
proposal rejects culture for. The operative statistic is therefore
min(S) ACROSS SPECIES, not mean(S) across the clade - this script reports
per-species minima and the fraction of genomes falling below a stated gate.

IMPORTANT - rRNA copy number and draft assemblies:
    Short-read assemblies collapse rRNA repeats, so counting rRNA-guide sites
    in a draft genome undercounts, often reporting 1-2 where the true copy
    number is 7-8. Run this with --closed-only to restrict to complete
    genomes when calibrating rRNA-guide copy number, and cross-check against
    rrnDB. The proposal already commits to this in §2.4; the flag enforces it.

Usage:
    python 24_count_sites.py --guides ../11_guides/candidates.fasta \\
        --genomes ../05_dereplication/output/dereplicated_genomes \\
        --out ../11_guides/site_counts.tsv

    # rRNA guides, complete genomes only, for copy-number calibration
    python 24_count_sites.py --guides ../11_guides/rrna_guides.fasta \\
        --closed-only --metadata ../01_metadata/Enterobacterales_high_quality_metadata.tsv
"""

import argparse
import re
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--guides", required=True, help="FASTA of 20 nt protospacers")
    ap.add_argument("--genomes", default="../05_dereplication/output/dereplicated_genomes")
    ap.add_argument("--metadata",
                    default="../02_genomes_raw/Clinical_Enterobacterales_discovery_8000.tsv")
    ap.add_argument("--out", default="../11_guides/site_counts.tsv")
    ap.add_argument("--gate", type=int, default=10,
                    help="Sites per genome assumed necessary to cross the FACS gate")
    ap.add_argument("--closed-only", action="store_true",
                    help="Restrict to complete genomes (rRNA copy-number calibration)")
    ap.add_argument("--threads", type=int, default=20)
    return ap.parse_args()


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
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


def load_genome(path):
    """Contigs joined by N runs so no match can straddle a contig boundary."""
    return ("N" * 30).join(seq for _, seq in read_fasta(path))


def count_with_pam(seq, motif):
    n, start, L = 0, 0, len(motif)
    while True:
        i = seq.find(motif, start)
        if i == -1:
            return n
        j = i + L
        if j + 3 <= len(seq) and seq[j + 1] == "G" and seq[j + 2] == "G":
            n += 1
        start = i + 1


def score_genome(path, guides):
    seq = load_genome(path)
    rev = revcomp(seq)
    name = Path(path).stem
    return name, [count_with_pam(seq, g) + count_with_pam(rev, g) for g in guides]


def main():
    args = parse_args()
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    guide_names, guides = [], []
    for h, s in read_fasta(args.guides):
        guide_names.append(h.split()[0])
        guides.append(s)
    print(f"Guides: {len(guides)}")

    genome_paths = sorted(Path(args.genomes).glob("*.fna"))
    print(f"Genomes: {len(genome_paths)}")

    # ------------------------------------------------------------------
    # Taxonomy and assembly level
    # ------------------------------------------------------------------
    tax, level = {}, {}
    meta_path = Path(args.metadata)
    if meta_path.exists():
        meta = pd.read_csv(meta_path, sep="\t", low_memory=False)
        for _, r in meta.iterrows():
            acc = r.get("ncbi_genbank_assembly_accession")
            if isinstance(acc, str):
                tax[acc] = str(r.get("gtdb_taxonomy", ""))
                level[acc] = str(r.get("ncbi_assembly_level", "")).lower()

    if args.closed_only:
        before = len(genome_paths)
        genome_paths = [p for p in genome_paths
                        if level.get(p.stem, "") in ("complete genome",)]
        print(f"  restricted to complete assemblies: {len(genome_paths)} of {before}")
        if not genome_paths:
            raise SystemExit(
                "No complete genomes found. Check ncbi_assembly_level is present "
                "in the metadata, or drop --closed-only.")

    # ------------------------------------------------------------------
    # Count
    # ------------------------------------------------------------------
    print(f"\nCounting sites on {args.threads} processes ...")
    worker = partial(score_genome, guides=guides)
    with Pool(args.threads) as pool:
        results = pool.map(worker, [str(p) for p in genome_paths], chunksize=8)

    df = pd.DataFrame([r[1] for r in results], columns=guide_names,
                      index=[r[0] for r in results])
    df.index.name = "Genome"
    df["S_total"] = df[guide_names].sum(axis=1)

    df["family"] = [re.search(r"f__([^;]*)", tax.get(g, "")).group(1)
                    if re.search(r"f__([^;]*)", tax.get(g, "")) else "unknown"
                    for g in df.index]
    df["species"] = [re.search(r"s__([^;]*)", tax.get(g, "")).group(1)
                     if re.search(r"s__([^;]*)", tax.get(g, "")) else "unknown"
                     for g in df.index]

    df.to_csv(out_path, sep="\t")

    # ------------------------------------------------------------------
    # Photon budget report
    # ------------------------------------------------------------------
    S = df["S_total"]
    print()
    print("PHOTON BUDGET  (sites per genome, whole panel)")
    print(f"  mean       {S.mean():.1f}")
    print(f"  median     {S.median():.0f}")
    print(f"  min        {S.min()}")
    print(f"  5th pct    {S.quantile(0.05):.0f}")
    print()
    print(f"  genomes with 0 sites:        {int((S == 0).sum())} "
          f"({100 * (S == 0).mean():.2f}%)  <- never recoverable")
    print(f"  genomes below gate ({args.gate}):     {int((S < args.gate).sum())} "
          f"({100 * (S < args.gate).mean():.2f}%)")

    # Per-species minimum is the gate-setting statistic (proposal §2.4)
    sp = df.groupby("species")["S_total"].agg(
        n="size", min="min", median="median", mean="mean")
    sp = sp[sp["n"] >= 3].sort_values("min")
    print()
    print("Weakest species by minimum sites (n>=3 genomes):")
    print(sp.head(20).round(1).to_string())

    fam = df.groupby("family")["S_total"].agg(
        n="size", min="min", median="median",
        pct_below_gate=lambda x: 100 * (x < args.gate).mean())
    print()
    print("By family:")
    print(fam.sort_values("median").round(1).to_string())

    print()
    print("Per-guide contribution (mean sites per genome):")
    contrib = df[guide_names].mean().sort_values(ascending=False)
    print(contrib.round(2).head(25).to_string())

    print()
    print("Saved:", out_path)
    print()
    print("Set the FACS gate against the weakest species above, not the mean.")
    print("Recovery r = fraction of genomes at or above gate; feed r into")
    print("19_select_panel.py, which uses it to compute the purity budget.")


if __name__ == "__main__":
    main()
