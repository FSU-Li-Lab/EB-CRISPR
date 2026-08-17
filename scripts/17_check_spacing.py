#!/usr/bin/env python3
"""
29_check_spacing.py

Map each guide onto its reference locus and flag sites too close together.

16_check_guide_overlap.py caught guides that SHARE sequence. This catches the
other half of the problem: two guides can have no sequence in common and still
compete, because a bound dCas9 protects roughly 23-30 bp of duplex - the R-loop
plus the protein footprint. Sites closer than that are effectively one site,
and counting both inflates the photon budget exactly as a redundant pair does.

This matters most for short loci. Three guides inside a 1.3 kb secY gene have
far less room than five spread across a 2.9 kb 23S.

The script also reports which reference sequences contain every guide, which is
a useful sanity check in itself: a guide that cannot be located in any single
reference is either mis-assembled in that genome or was called from a
divergent copy.

Usage:
    python3 29_check_spacing.py \\
        --guides ../11_guides/panel_nonredundant.fasta \\
        --locus 23S ../10_rrna/23S_one_per_genome.fa \\
        --locus secY ../10_core_genes/secY.fna
"""

import argparse

COMP = str.maketrans("ACGTN", "TGCAN")


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--guides", required=True)
    ap.add_argument("--locus", nargs=2, action="append", metavar=("PREFIX", "FASTA"),
                    required=True,
                    help="Guide-name prefix and the FASTA of that locus")
    ap.add_argument("--min-spacing", type=int, default=30,
                    help="Minimum bp between adjacent sites (dCas9 footprint)")
    ap.add_argument("--max-refs", type=int, default=2000,
                    help="Reference sequences to scan when picking the best one")
    return ap.parse_args()


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path, limit=None):
    name, seq, out = None, [], {}
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name:
                out[name] = "".join(seq).upper()
                if limit and len(out) >= limit:
                    return out
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line)
    if name:
        out[name] = "".join(seq).upper()
    return out


def locate(ref, guide):
    """Return (position, strand) of guide in ref, or None."""
    i = ref.find(guide)
    if i >= 0:
        return i, "+"
    i = ref.find(revcomp(guide))
    if i >= 0:
        return i, "-"
    return None


def main():
    args = parse_args()
    guides = read_fasta(args.guides)
    print(f"Guides: {len(guides)}\n")

    any_warning = False

    for prefix, path in args.locus:
        subset = {n: s for n, s in guides.items() if n.startswith(prefix)}
        if not subset:
            print(f"{prefix}: no guides with this prefix, skipping\n")
            continue

        refs = read_fasta(path, limit=args.max_refs)
        if not refs:
            print(f"{prefix}: could not read {path}\n")
            continue

        # Pick the reference containing the most guides - a genome missing one
        # would otherwise give a misleading spacing picture.
        best_ref, best_hits = None, -1
        for rname, rseq in refs.items():
            hits = sum(1 for s in subset.values() if locate(rseq, s))
            if hits > best_hits:
                best_ref, best_hits = rname, hits
            if hits == len(subset):
                break

        rseq = refs[best_ref]
        print(f"=== {prefix} ===")
        print(f"reference: {best_ref}  ({len(rseq)} bp)")
        print(f"guides located: {best_hits}/{len(subset)}\n")

        located = []
        for n, s in sorted(subset.items()):
            hit = locate(rseq, s)
            if hit:
                located.append((hit[0], n, hit[1]))
            else:
                print(f"  NOT FOUND in this reference: {n}")

        print(f"  {'pos':>7}  {'guide':<16} strand")
        for pos, n, strand in sorted(located):
            print(f"  {pos:>7}  {n:<16} {strand}")

        ordered = sorted(located)
        print()
        for (p1, n1, _), (p2, n2, _) in zip(ordered, ordered[1:]):
            gap = p2 - p1
            if gap < args.min_spacing:
                any_warning = True
                print(f"  WARNING  {n1} and {n2} are {gap} bp apart "
                      f"(< {args.min_spacing}) - sterically competing, "
                      f"count as ONE site")
            else:
                print(f"  ok       {n1} -> {n2}: {gap} bp")
        print()

    if any_warning:
        print("Remove one guide from each competing pair and recount the photon")
        print("budget with 18_count_sites.py.")
    else:
        print("All sites adequately separated - signal should be additive.")


if __name__ == "__main__":
    main()
