#!/usr/bin/env python3
"""
26_summarise_offtargets.py

Turn raw cas-offinder output into per-guide verdicts and the summary table
19_select_panel.py consumes.

Raw hit counts are close to meaningless on their own. 300,000 hits could be one
promiscuous guide or all of them, and it matters enormously which. What
determines whether a guide is usable is:

    - does any background genome carry enough sites to CLEAR THE GATE
    - how abundant is the taxon that carries them

A GATE CAVEAT SPECIFIC TO rRNA GUIDES
-------------------------------------
Earlier reasoning held that the gate protects you, because targets carry 7-8
rRNA sites while a background off-target is typically single-copy and sits
below threshold. That argument does NOT hold for rRNA guides. Background
bacteria have multi-copy rRNA operons too - a Vibrio or Haemophilus genome
carries roughly as many rRNA operons as an Enterobacterales one. An rRNA guide
that cross-reacts therefore cross-reacts 7 times over and clears the gate just
as easily as a true target.

So the copy-number margin only buys specificity for SINGLE-COPY guides against
single-copy off-targets. rRNA guides must be clean outright - near-zero
tolerance, judged on presence rather than on clearing a threshold. This script
reports both views: --gate filtering for core-gene guides, and raw presence for
rRNA guides.

Usage:
    python3 26_summarise_offtargets.py \\
        --offtarget-dir ../11_guides/offtarget_tier1 \\
        --guides ../11_guides/all_candidates.fasta \\
        --outdir ../11_guides/offtarget_tier1
"""

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--offtarget-dir", required=True)
    ap.add_argument("--guides", required=True)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--gate", type=int, default=10,
                    help="Sites in one background genome needed to clear the gate")
    ap.add_argument("--locus-map", default=None,
                    help="Optional TSV: guide name, locus. Otherwise inferred "
                         "from the candidate FASTA filenames if names carry them.")
    return ap.parse_args()


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


def load_hits(path):
    """cas-offinder output, tolerant of header lines and column drift.
    Returns list of (query, contig)."""
    out = []
    if not Path(path).exists():
        return out
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 2:
                continue
            # Some builds prepend an id column; the query is the first field
            # made only of ACGTNRYacgtn and at least 20 chars.
            q, c = f[0], f[1]
            if not set(q.upper()) <= set("ACGTNRY"):
                if len(f) >= 3:
                    q, c = f[1], f[2]
                else:
                    continue
            out.append((q.upper(), c))
    return out


def main():
    args = parse_args()
    od = Path(args.offtarget_dir)
    outdir = Path(args.outdir or args.offtarget_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    guides = {h.split()[0]: s for h, s in read_fasta(args.guides)}
    seq_to_names = defaultdict(list)
    seed_to_names = defaultdict(list)
    for n, s in guides.items():
        seq_to_names[s].append(n)
        seed_to_names[s[-12:]].append(n)
    print(f"Guides: {len(guides)}")

    # ------------------------------------------------------------------
    # Parse both passes
    # ------------------------------------------------------------------
    # full: query is {protospacer}NNN   -> strip trailing N
    # seed: query is NNNNNNNN{seed}NNN  -> strip both ends
    per_guide_genome = defaultdict(lambda: defaultdict(int))   # guide -> acc -> sites
    genome_genus = {}
    seed_hit_guides = set()

    strict = load_hits(od / "offtargets_strict.txt")
    print(f"Full-length hit rows: {len(strict)}")
    for q, contig in strict:
        proto = q.rstrip("N")
        names = seq_to_names.get(proto, [])
        if not names:
            continue
        parts = contig.split("|")
        genus = parts[0] if len(parts) >= 3 else "unknown"
        acc = parts[1] if len(parts) >= 3 else contig
        genome_genus[acc] = genus
        for n in names:
            per_guide_genome[n][acc] += 1

    seed = load_hits(od / "offtargets_seed.txt")
    print(f"Seed hit rows:        {len(seed)}")
    for q, contig in seed:
        s = q.strip("N")
        for n in seed_to_names.get(s, []):
            seed_hit_guides.add(n)

    # ------------------------------------------------------------------
    # Per-guide verdict
    # ------------------------------------------------------------------
    rows = []
    for name, proto in guides.items():
        hits = per_guide_genome.get(name, {})
        n_genomes = len(hits)
        over_gate = {a: c for a, c in hits.items() if c >= args.gate}
        genera = {genome_genus.get(a, "unknown") for a in hits}
        rows.append({
            "Guide": name,
            "Protospacer": proto,
            "bg_genomes_hit": n_genomes,
            "bg_genomes_over_gate": len(over_gate),
            "bg_genera_hit": len(genera),
            "total_sites": sum(hits.values()),
            "seed_hit": name in seed_hit_guides,
            "top_genera": ";".join(sorted(genera)[:5]),
        })

    df = pd.DataFrame(rows).sort_values(
        ["bg_genomes_hit", "total_sites"])
    df.to_csv(outdir / "per_guide_offtargets.tsv", sep="\t", index=False)

    clean_full = df[df["bg_genomes_hit"] == 0]
    clean_seed = df[~df["seed_hit"]]
    clean_both = df[(df["bg_genomes_hit"] == 0) & (~df["seed_hit"])]
    under_gate = df[(df["bg_genomes_over_gate"] == 0) & (~df["seed_hit"])]

    print()
    print("VERDICT")
    print(f"  guides with zero full-length background hits: {len(clean_full)}")
    print(f"  guides with zero seed hits:                   {len(clean_seed)}")
    print(f"  clean on both:                                {len(clean_both)}")
    print(f"  no background genome over gate ({args.gate}), no seed hit: {len(under_gate)}")
    print()
    print("  The last row is the usable set for SINGLE-COPY guides only.")
    print("  For rRNA guides use 'clean on both' - background rRNA is multi-copy,")
    print("  so a cross-reacting rRNA guide clears the gate as easily as a target.")

    if len(clean_both):
        with open(outdir / "surviving_candidates.fasta", "w") as fh:
            for r in clean_both.itertuples():
                fh.write(f">{r.Guide}\n{r.Protospacer}\n")
        print(f"\n  wrote surviving_candidates.fasta ({len(clean_both)})")

    if len(under_gate):
        with open(outdir / "surviving_under_gate.fasta", "w") as fh:
            for r in under_gate.itertuples():
                fh.write(f">{r.Guide}\n{r.Protospacer}\n")
        print(f"  wrote surviving_under_gate.fasta ({len(under_gate)})")

    print()
    print("Worst offenders:")
    print(df.sort_values("bg_genomes_hit", ascending=False)
            [["Guide", "bg_genomes_hit", "bg_genera_hit", "total_sites", "seed_hit"]]
            .head(10).to_string(index=False))

    print()
    print("Best candidates:")
    print(df[["Guide", "bg_genomes_hit", "bg_genomes_over_gate",
              "bg_genera_hit", "seed_hit"]].head(15).to_string(index=False))

    # ------------------------------------------------------------------
    # Input for 19_select_panel.py
    # ------------------------------------------------------------------
    summary = []
    for name, hits in per_guide_genome.items():
        by_genus = defaultdict(int)
        for acc, c in hits.items():
            if c >= args.gate:
                by_genus[genome_genus.get(acc, "unknown")] = max(
                    by_genus[genome_genus.get(acc, "unknown")], c)
        for genus, c in by_genus.items():
            summary.append({"Guide_seq": guides[name], "Taxon": genus, "N_sites": c})

    sdf = pd.DataFrame(summary, columns=["Guide_seq", "Taxon", "N_sites"])
    sdf.to_csv(outdir / "summary.tsv", sep="\t", index=False)
    print(f"\nSaved: {outdir / 'summary.tsv'}  ({len(sdf)} guide-taxon pairs)")
    print("Saved:", outdir / "per_guide_offtargets.tsv")


if __name__ == "__main__":
    main()
