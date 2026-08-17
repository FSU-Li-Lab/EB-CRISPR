#!/usr/bin/env python3
"""
25_select_panel.py

Turn milestone M1.3 (<=5% non-EB carryover) into a guide-design criterion, then
assemble the panel that meets it.

THE PURITY BUDGET
-----------------
The proposal states that M1.3 "is set at guide design rather than at the
sorter". That is right, and it has an arithmetic consequence worth making
explicit, because it is much sharper than it looks.

Let
    a  = Enterobacterales fraction of the community          (~0.01)
    r  = fraction of EB cells labeled above the gate         (recovery)
    e  = fraction of non-EB cells labeled above the gate     (false label rate)
    P  = required purity of the sorted population            (0.95)

    purity = a*r / (a*r + (1-a)*e)  >=  P

    =>   e  <=  a * r * (1 - P) / (P * (1 - a))

At a = 0.01, r = 0.5, P = 0.95 this gives e <= 2.7e-4. Fewer than about one in
3,700 non-target cells may be labeled. That is the real specificity
specification, and it is roughly three orders of magnitude stricter than
"no off-targets in the eight dominant gut genera".

Two things follow.

1. The counter-screen must be ABUNDANCE-WEIGHTED. An off-target in
   Bacteroides at 20% of the community consumes the entire budget on its own,
   even at 1% labeling efficiency. The same off-target in a genus at 0.01%
   abundance is irrelevant. Screening a flat list of genomes treats these as
   equivalent; they differ by three orders of magnitude. Weight by the
   abundance profile of the samples you will actually run - once Task 4 bulk
   metagenomics on the twenty FSU specimens exists, use that profile rather
   than any published average.

2. Recovery and purity trade against each other, and the trade is favourable.
   Because e scales linearly with r, raising the gate to exclude weakly
   labeled non-target cells costs recovery but buys purity proportionally.
   The script reports the frontier so the gate can be placed deliberately.

PANEL ASSEMBLY
--------------
Guides are selected greedily to maximise r - the fraction of target genomes at
or above the gate - subject to the summed off-target budget. Maximising r
rather than mean signal is the point: a guide that adds brightness to already
bright genomes adds nothing, while one that rescues the weakest species raises
recovery directly. This is also what protects M1.2, since a species sitting
below the gate is partially recovered and partial recovery is compositional
distortion.

Note on panel size: off-target surface grows linearly with the number of
guides, and so does signal, so panel size alone does not improve the
signal-to-background ratio. It improves the WORST CASE. Add guides to lift
weak species, not to lift the mean.

Usage:
    python 25_select_panel.py \\
        --site-counts ../11_guides/site_counts.tsv \\
        --guides ../11_guides/candidates.fasta \\
        --offtargets ../11_guides/offtarget/summary.tsv \\
        --abundances ../09_background/gut_abundance.tsv

Build the offtarget summary from cas-offinder output (guide seq, taxon, sites):
    awk 'BEGIN{OFS="\\t"} {print $1, $2, 1}' offtargets_seed.txt \\
      | sort | uniq -c | awk 'BEGIN{OFS="\\t"}{print $2,$3,$1}' \\
      > summary.tsv
    (then map contig -> taxon for column 2)
"""

import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--site-counts", default="../11_guides/site_counts.tsv")
    ap.add_argument("--guides", required=True)
    ap.add_argument("--offtargets", default=None,
                    help="TSV: Guide_seq, Taxon, N_sites")
    ap.add_argument("--abundances", default=None,
                    help="TSV: Taxon, Relative_abundance (fractions summing to ~1)")
    ap.add_argument("--out", default="../11_guides/panel.tsv")
    ap.add_argument("--eb-abundance", type=float, default=0.01)
    ap.add_argument("--purity", type=float, default=0.95)
    ap.add_argument("--gate", type=int, default=10)
    ap.add_argument("--max-guides", type=int, default=20)
    ap.add_argument("--label-efficiency", type=float, default=1.0,
                    help="P(non-target cell labeled | it carries an off-target site). "
                         "1.0 is the conservative default.")
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


def purity_budget(a, r, P):
    """Max tolerable false-label rate among non-target cells."""
    return a * r * (1 - P) / (P * (1 - a))


def main():
    args = parse_args()

    df = pd.read_csv(args.site_counts, sep="\t", index_col="Genome")
    meta_cols = [c for c in ("S_total", "family", "species") if c in df.columns]
    guide_cols = [c for c in df.columns if c not in meta_cols]
    print(f"Genomes: {len(df)}   candidate guides: {len(guide_cols)}")

    name_to_seq = {h.split()[0]: s for h, s in read_fasta(args.guides)}

    # ------------------------------------------------------------------
    # Off-target load per guide, abundance weighted
    # ------------------------------------------------------------------
    e_guide = {g: 0.0 for g in guide_cols}

    if args.offtargets and Path(args.offtargets).exists():
        ot = pd.read_csv(args.offtargets, sep="\t",
                         names=["Guide_seq", "Taxon", "N_sites"], header=0)

        if args.abundances and Path(args.abundances).exists():
            ab = pd.read_csv(args.abundances, sep="\t",
                             names=["Taxon", "Relative_abundance"], header=0)
            weights = dict(zip(ab["Taxon"], ab["Relative_abundance"]))
            total_w = sum(weights.values())
            print(f"Abundance profile: {len(weights)} taxa, sums to {total_w:.3f}")
        else:
            taxa = ot["Taxon"].unique()
            weights = {t: 1.0 / len(taxa) for t in taxa}
            print()
            print("WARNING: no abundance profile supplied. Using uniform weights.")
            print("  This UNDERSTATES risk from dominant taxa by orders of magnitude")
            print("  and the resulting budget should not be trusted. Supply a real")
            print("  profile - ideally from bulk metagenomics of your own samples.")

        seq_to_name = {v: k for k, v in name_to_seq.items()}
        for _, row in ot.iterrows():
            gname = seq_to_name.get(row["Guide_seq"])
            if gname in e_guide:
                e_guide[gname] += (weights.get(row["Taxon"], 0.0)
                                   * args.label_efficiency)
    else:
        print()
        print("WARNING: no off-target file. Selecting on signal only.")
        print("  The resulting panel is NOT specificity-validated.")

    # ------------------------------------------------------------------
    # Purity budget
    # ------------------------------------------------------------------
    a, P = args.eb_abundance, args.purity
    print()
    print("PURITY BUDGET")
    print(f"  EB abundance a = {a}, required purity P = {P}")
    print()
    print(f"  {'recovery r':>12} {'max e':>12} {'1 in N':>12}")
    for r in (0.2, 0.3, 0.5, 0.7, 0.9):
        e = purity_budget(a, r, P)
        print(f"  {r:>12.2f} {e:>12.2e} {1 / e:>12,.0f}")

    # ------------------------------------------------------------------
    # Greedy panel assembly
    # ------------------------------------------------------------------
    e_max = purity_budget(a, 0.5, P)      # provisional, refined below
    selected, e_used = [], 0.0
    current = pd.Series(0, index=df.index, dtype=int)

    print()
    print(f"Assembling panel (provisional budget e <= {e_max:.2e}) ...")

    # Objective is (recovery, progress), compared lexicographically.
    # "progress" is mean(min(S, gate)): how far genomes have moved TOWARD the
    # gate, saturating once they clear it. Without it the greedy is myopic. A
    # guide lifting a weak family from 1 site to 3 changes recovery by zero at
    # a gate of 5, so it would never be selected - even though it is precisely
    # what a second guide needs in order to push that family over.
    def objective(S):
        return (float((S >= args.gate).mean()),
                float(S.clip(upper=args.gate).mean()))

    for step in range(args.max_guides):
        best, best_obj = None, (-1.0, -1.0)
        for g in guide_cols:
            if g in selected:
                continue
            if e_used + e_guide[g] > e_max:
                continue
            obj = objective(current + df[g])
            if obj > best_obj or (obj == best_obj and best is not None
                                  and e_guide[g] < e_guide[best]):
                best, best_obj = g, obj

        if best is None:
            print("  budget exhausted")
            break

        if step > 0 and best_obj <= objective(current):
            print("  no further gain in recovery or progress toward gate")
            break

        best_r = best_obj[0]

        selected.append(best)
        e_used += e_guide[best]
        current = current + df[best]
        print(f"  + {best:<12} r = {best_r:.4f}   sum(e) = {e_used:.2e}")

    if not selected:
        raise SystemExit("No guides selected - budget too tight or no candidates.")

    # ------------------------------------------------------------------
    # Final evaluation, with the budget recomputed at the achieved r
    # ------------------------------------------------------------------
    r_final = float((current >= args.gate).mean())
    e_max_final = purity_budget(a, r_final, P)
    purity_est = (a * r_final) / (a * r_final + (1 - a) * e_used) if e_used > 0 else 1.0

    print()
    print("PANEL")
    for g in selected:
        print(f"  {g:<12} {name_to_seq.get(g, '?')}  "
              f"mean sites/genome {df[g].mean():.2f}  e = {e_guide[g]:.2e}")

    print()
    print(f"  recovery r              {r_final:.4f}")
    print(f"  summed off-target e     {e_used:.2e}")
    print(f"  budget at this r        {e_max_final:.2e}"
          f"   {'OK' if e_used <= e_max_final else 'EXCEEDED'}")
    print(f"  predicted purity        {100 * purity_est:.2f}%   (M1.3 needs >=95%)")
    print()
    print(f"  sites per genome S:  median {current.median():.0f}, "
          f"min {current.min()}, 5th pct {current.quantile(0.05):.0f}")
    print(f"  -> this is the [S] bracket in proposal section 2.5")

    if "family" in df.columns:
        fam = pd.DataFrame({"S": current, "family": df["family"]})
        summary = fam.groupby("family")["S"].agg(
            n="size", median="median", min="min",
            recovery=lambda x: (x >= args.gate).mean())
        print()
        print("Recovery by family:")
        print(summary.sort_values("recovery").round(3).to_string())

        weak = summary[summary["recovery"] < 0.9]
        if len(weak):
            print()
            print("Families below 90% recovery. Options, in the order the proposal")
            print("already anticipates them:")
            print("  1. add family-specific guides (raises e, check budget)")
            print("  2. lower the gate (raises r AND raises e proportionally)")
            print("  3. restrict the claim to Enterobacteriaceae, per the risk table")

    out = pd.DataFrame({
        "Guide": selected,
        "Protospacer": [name_to_seq.get(g, "") for g in selected],
        "Mean_sites_per_genome": [df[g].mean() for g in selected],
        "Offtarget_e": [e_guide[g] for g in selected],
    })
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print()
    print("Saved:", args.out)


if __name__ == "__main__":
    main()
