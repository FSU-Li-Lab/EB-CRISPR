#!/usr/bin/env python3
"""
20_build_negative_set.py

Build the background (negative) genome set that defines guide specificity.

You already have bac120_metadata.tsv on disk, so most of this costs nothing
extra to assemble. The set is tiered because the tiers do different jobs:

  Tier 1  HARD NEGATIVES - Gammaproteobacteria outside Enterobacterales
          (Pasteurellales, Pseudomonadales, Vibrionales, Aeromonadales,
          Moraxellales, Xanthomonadales, Burkholderiales*). These are where
          false positives will actually come from. Nucleotide identity to
          Enterobacterales in rRNA and housekeeping loci is high enough that
          a guide can easily cross-react. Screen against every one of these.

  Tier 2  GUT BACKGROUND - the dominant phyla in stool. Large, phylogenetically
          distant, unlikely to cross-react, but they are the bulk of the DNA
          in your sample so they must be screened.

  Tier 3  BROAD - remaining bacterial diversity. Optional; use if you want to
          claim specificity beyond the gut context.

Not covered here, add separately:
  - UHGG / MGnify human-gut catalogue (species reps), for gut strains missing
    from GTDB:
    ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/
  - GRCh38 and T2T-CHM13, if host DNA survives your extraction
  - Archaea (Methanobrevibacter smithii is common in gut) - set --archaea

Usage:
    python 20_build_negative_set.py
    python 20_build_negative_set.py --tiers 1 2 --max-per-tier 3000
"""

import argparse
from pathlib import Path

import pandas as pd

# Orders within Gammaproteobacteria that are NOT Enterobacterales.
# Left as a class-level match rather than an explicit list, because GTDB
# reshuffles Gammaproteobacterial orders between releases and a hardcoded
# list silently goes stale.
HARD_CLASS = "c__Gammaproteobacteria"
TARGET_ORDER = "o__Enterobacterales"

# Gammaproteobacterial orders adjacent to Enterobacterales. These are where
# cross-reactivity actually originates: close enough that rRNA and housekeeping
# loci share long conserved stretches. Never subsample these.
#
# Note that GTDB folds the old Betaproteobacteria (Burkholderiales and
# relatives) into c__Gammaproteobacteria. Those are environmental organisms,
# numerous in GTDB and largely absent from gut samples, so a random subsample
# of "Gammaproteobacteria" is dominated by taxa that pose almost no risk while
# diluting the ones that do.
NEAR_ORDERS = [
    "o__Pasteurellales",
    "o__Vibrionales",
    "o__Aeromonadales",
    "o__Alteromonadales",
    "o__Oceanospirillales",
    "o__Cardiobacteriales",
    "o__Orbales",
    "o__Enterobacterales_A",
]

# Orders that dominate the human gut. Risk scales with abundance, so these are
# kept in full regardless of how many genomes that means.
GUT_ORDERS = [
    "o__Bacteroidales",
    "o__Lachnospirales",
    "o__Oscillospirales",
    "o__Christensenellales",
    "o__Erysipelotrichales",
    "o__Veillonellales",
    "o__Selenomonadales",
    "o__Peptostreptococcales",
    "o__Clostridiales",
    "o__Bifidobacteriales",
    "o__Coriobacteriales",
    "o__Actinomycetales",
    "o__Lactobacillales",
    "o__Campylobacterales",
    "o__Verrucomicrobiales",
    "o__Desulfovibrionales",
    "o__Fusobacteriales",
    "o__Enterobacterales_B",
]

GUT_PHYLA = [
    "p__Bacteroidota",
    "p__Bacillota",        # GTDB name for Firmicutes
    "p__Bacillota_A",
    "p__Bacillota_B",
    "p__Bacillota_C",
    "p__Firmicutes",       # older GTDB releases
    "p__Firmicutes_A",
    "p__Firmicutes_B",
    "p__Firmicutes_C",
    "p__Actinobacteriota",
    "p__Verrucomicrobiota",
    "p__Desulfobacterota",
    "p__Desulfobacterota_A",
    "p__Fusobacteriota",
    "p__Campylobacterota",
    "p__Synergistota",
    "p__Spirochaetota",
]


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", default="../01_metadata/bac120_metadata.tsv")
    ap.add_argument("--outdir", default="../09_background")
    ap.add_argument("--tiers", nargs="+", type=int, default=[1, 2],
                    choices=[1, 2, 3])
    ap.add_argument("--cap-per-order", type=int, default=150,
                    help="Max genomes per GTDB order. Orders in NEAR_ORDERS "
                         "and GUT_ORDERS are exempt and always kept in full.")
    ap.add_argument("--max-per-tier", type=int, default=None,
                    help="Hard ceiling per tier, applied after the per-order "
                         "cap. Protected orders are never trimmed by it.")
    ap.add_argument("--min-completeness", type=float, default=90.0)
    ap.add_argument("--max-contamination", type=float, default=5.0)
    ap.add_argument("--seed", type=int, default=42)
    return ap.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Reading {args.metadata} ...")
    df = pd.read_csv(args.metadata, sep="\t", low_memory=False)
    print(f"  total GTDB genomes: {len(df)}")

    # ------------------------------------------------------------------
    # Species representatives only. The background does not need strain
    # depth - it needs breadth. One genome per species keeps the
    # off-target search tractable.
    # ------------------------------------------------------------------
    if "gtdb_representative" in df.columns:
        reps = df[df["gtdb_representative"].astype(str).str.lower().isin(["t", "true"])]
    else:
        reps = df[df["gtdb_genome_representative"] == df["accession"]]
    print(f"  species representatives: {len(reps)}")

    reps = reps[
        (reps["checkm2_completeness"] >= args.min_completeness) &
        (reps["checkm2_contamination"] <= args.max_contamination)
    ].copy()
    print(f"  after quality filter: {len(reps)}")

    tax = reps["gtdb_taxonomy"].fillna("")

    # Everything Enterobacterales is excluded from the background by
    # definition - these are the positives.
    is_target = tax.str.contains(TARGET_ORDER, na=False)
    print(f"  excluded as target (Enterobacterales): {int(is_target.sum())}")
    bg = reps[~is_target].copy()
    bgtax = bg["gtdb_taxonomy"].fillna("")

    # ------------------------------------------------------------------
    # Tier assignment
    # ------------------------------------------------------------------
    tier1 = bg[bgtax.str.contains(HARD_CLASS, na=False)]

    gut_pattern = "|".join(GUT_PHYLA)
    tier2 = bg[bgtax.str.contains(gut_pattern, na=False)]

    assigned = set(tier1["accession"]) | set(tier2["accession"])
    tier3 = bg[~bg["accession"].isin(assigned)]

    tiers = {1: tier1, 2: tier2, 3: tier3}
    names = {1: "hard_negatives_gammaproteobacteria",
             2: "gut_background",
             3: "broad_bacteria"}

    print()
    selected = []
    keep_all = set(NEAR_ORDERS) | set(GUT_ORDERS)

    for t in sorted(tiers):
        sub = tiers[t]
        mark = "*" if t in args.tiers else " "
        print(f" {mark} tier {t} ({names[t]}): {len(sub)} genomes")
        if t not in args.tiers:
            continue

        sub = sub.copy()
        sub["order"] = sub["gtdb_taxonomy"].str.extract(r"(o__[^;]+)")[0].fillna("unknown")

        # Cap PER ORDER rather than subsampling the tier at random. A global
        # random draw fills its quota proportionally, so the largest orders in
        # GTDB crowd out the ones that matter. Risk here is driven by
        # phylogenetic proximity and by gut abundance, neither of which
        # correlates with how many genomes an order happens to have deposited.
        kept = []
        for order, grp in sub.groupby("order"):
            if order in keep_all or args.cap_per_order is None:
                kept.append(grp)
            elif len(grp) > args.cap_per_order:
                kept.append(grp.sample(n=args.cap_per_order, random_state=args.seed))
            else:
                kept.append(grp)
        sub = pd.concat(kept, ignore_index=True)

        if args.max_per_tier and len(sub) > args.max_per_tier:
            protected = sub[sub["order"].isin(keep_all)]
            rest = sub[~sub["order"].isin(keep_all)]
            room = max(args.max_per_tier - len(protected), 0)
            if room and len(rest) > room:
                rest = rest.sample(n=room, random_state=args.seed)
            sub = pd.concat([protected, rest], ignore_index=True)
            print(f"      capped to {len(sub)} "
                  f"({len(protected)} protected, {len(rest)} sampled)")
        else:
            print(f"      after per-order cap: {len(sub)}")

        sub["background_tier"] = t
        selected.append(sub)

    if not selected:
        raise SystemExit("No tiers selected")

    combined = pd.concat(selected, ignore_index=True)

    # ------------------------------------------------------------------
    # Report composition - this is what you cite as your specificity claim
    # ------------------------------------------------------------------
    order = combined["gtdb_taxonomy"].str.extract(r"(o__[^;]+)")[0]
    print()
    print("Top orders in background:")
    print(order.value_counts().head(15).to_string())

    phylum = combined["gtdb_taxonomy"].str.extract(r"(p__[^;]+)")[0]
    print()
    print("Phyla in background:")
    print(phylum.value_counts().head(15).to_string())

    # ------------------------------------------------------------------
    # Write
    # ------------------------------------------------------------------
    meta_out = outdir / "background_metadata.tsv"
    combined.to_csv(meta_out, sep="\t", index=False)

    acc = combined["ncbi_genbank_assembly_accession"].dropna().drop_duplicates()
    acc_out = outdir / "background_accessions.txt"
    acc.to_csv(acc_out, index=False, header=False)

    # Per-tier lists, so you can screen hard negatives first and cheaply
    for t in args.tiers:
        sub = combined[combined["background_tier"] == t]
        a = sub["ncbi_genbank_assembly_accession"].dropna().drop_duplicates()
        a.to_csv(outdir / f"background_tier{t}_accessions.txt",
                 index=False, header=False)

    print()
    print(f"Background genomes: {len(combined)}  ({len(acc)} with NCBI accessions)")
    print("Saved:")
    print(" ", meta_out)
    print(" ", acc_out)
    print()
    print("Download with:")
    print(f"  datasets download genome accession --inputfile {acc_out} \\")
    print(f"      --include genome --filename {outdir}/background.zip")


if __name__ == "__main__":
    main()
