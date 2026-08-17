#!/usr/bin/env python3
"""
27_check_guide_overlap.py

Find guides that target the same locus and therefore cannot bind at once.

dCas9 unwinds the DNA duplex to form an R-loop, so an occupied target site is
unavailable to any other complex - including one whose protospacer lies on the
opposite strand. Two guides pointing at the same locus from opposite strands
are therefore ONE binding site, not two, and counting both inflates the photon
budget.

This is easy to miss because such a pair looks completely dissimilar when
compared directly: the relationship only appears after reverse-complementing
one of them. Their tell in the coverage tables is identical presence/absence
across every genome.

Overlap is reported as the longest shared substring between guide A (or its
reverse complement) and guide B. Anything above ~10 nt means the two sit at
essentially the same locus; below that they may still be close enough to
interfere sterically, which position mapping rather than sequence comparison
will settle.

Usage:
    python3 27_check_guide_overlap.py --guides ../11_guides/clean10.fasta
    python3 27_check_guide_overlap.py --guides ../11_guides/clean10.fasta \\
        --out ../11_guides/panel_nonredundant.fasta
"""

import argparse
from itertools import combinations

COMP = str.maketrans("ACGTN", "TGCAN")


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--guides", required=True)
    ap.add_argument("--min-overlap", type=int, default=10)
    ap.add_argument("--out", default=None,
                    help="Write a de-duplicated FASTA keeping one guide per site")
    ap.add_argument("--site-counts", default=None,
                    help="Optional site_counts TSV; used to keep the guide with "
                         "broader genome coverage from each redundant pair")
    return ap.parse_args()


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    name, seq, out = None, [], {}
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name:
                out[name] = "".join(seq).upper()
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line)
    if name:
        out[name] = "".join(seq).upper()
    return out


def longest_common_substring(a, b):
    best = 0
    for i in range(len(a)):
        for j in range(i + best + 1, len(a) + 1):
            if a[i:j] in b:
                best = j - i
            else:
                break
    return best


def main():
    args = parse_args()
    guides = read_fasta(args.guides)
    print(f"Guides: {len(guides)}\n")

    coverage = {}
    if args.site_counts:
        import pandas as pd
        d = pd.read_csv(args.site_counts, sep="\t")
        for g in guides:
            if g in d.columns:
                coverage[g] = float((d[g] > 0).mean())

    conflicts = []
    for a, b in combinations(sorted(guides), 2):
        sa, sb = guides[a], guides[b]
        fwd = longest_common_substring(sa, sb)
        rev = longest_common_substring(revcomp(sa), sb)
        best = max(fwd, rev)
        if best >= args.min_overlap:
            conflicts.append((a, b, best, "same strand" if fwd >= rev
                              else "opposite strands"))

    if not conflicts:
        print("No competing guides found - all target distinct sites.")
    else:
        print(f"{'guide A':<16} {'guide B':<16} {'overlap':>8}  orientation")
        for a, b, n, orient in sorted(conflicts, key=lambda x: -x[2]):
            print(f"{a:<16} {b:<16} {n:>6} nt  {orient}")
        print()
        print("Guides sharing a locus bind mutually exclusively. Count ONE of")
        print("each pair toward the photon budget.")

    # ------------------------------------------------------------------
    # De-duplicate: union-find over conflicts, keep best per group
    # ------------------------------------------------------------------
    parent = {g: g for g in guides}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b, _, _ in conflicts:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    groups = {}
    for g in guides:
        groups.setdefault(find(g), []).append(g)

    keep = []
    for root, members in sorted(groups.items()):
        if len(members) == 1:
            keep.append(members[0])
        else:
            best = max(members, key=lambda m: (coverage.get(m, 0), m))
            keep.append(best)
            dropped = [m for m in members if m != best]
            print(f"\n  site group {sorted(members)} -> keeping {best}"
                  + (f" (coverage {coverage[best]:.2f})" if best in coverage else "")
                  + f", dropping {dropped}")

    print(f"\nDistinct binding sites: {len(keep)} of {len(guides)} guides")

    if args.out:
        with open(args.out, "w") as fh:
            for g in sorted(keep):
                fh.write(f">{g}\n{guides[g]}\n")
        print(f"Saved: {args.out}")
        print("\nRecount the photon budget with 18_count_sites.py on this file -")
        print("the previous S was inflated by the redundant pairs.")


if __name__ == "__main__":
    main()
