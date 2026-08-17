#!/usr/bin/env python3
"""
22_scan_conserved_grna.py

Find 20 nt protospacers with an adjacent NGG PAM that are conserved across the
target set. This is the actual guide-discovery step.

Design model for a dCas9 BINDING assay (not cleavage):

  dCas9 binding and Cas9 cleavage have different specificity profiles.
  Cleavage needs near-perfect complementarity. Binding is more permissive at
  the PAM-distal end and is dominated by the PAM plus the PAM-proximal seed
  (~10-12 nt). For a biosensor this cuts both ways:

    - It HELPS coverage: a guide can still bind an Enterobacterales genome
      carrying 1-2 PAM-distal substitutions, so your effective breadth across
      the order is wider than exact matching suggests.

    - It HURTS specificity: a background organism needs less similarity to
      produce a signal than it would to be cleaved. Off-target screening must
      therefore be run with a GENEROUS mismatch budget (see 23_offtarget).

  So this script reports two numbers per candidate:

    cov_exact  fraction of targets with a perfect 20 nt + NGG match
    cov_seed   fraction with a perfect PAM + 12 nt seed match, PAM-distal
               end unconstrained. This is the optimistic upper bound on
               binding coverage.

  Use cov_exact to rank and cov_seed to understand how much of the gap a
  degenerate or multiplexed guide set could close. Do not report cov_seed as
  though it were validated coverage - binding affinity falls off with
  PAM-distal mismatches, it does not stay flat.

Candidates are enumerated from one reference sequence per family, not just
from E. coli, so that windows absent from Enterobacteriaceae are still
considered.

Usage:
    # after: mafft --auto 16S_one_per_genome.fa > 16S_aln.fa  (alignment
    # optional - this script works on unaligned sequence)
    python 22_scan_conserved_grna.py \\
        --fasta ../10_rrna/16S_one_per_genome.fa \\
        --metadata ../02_genomes_raw/Clinical_Enterobacterales_discovery_8000.tsv \\
        --out ../11_guides/16S_candidates.tsv
"""

import argparse
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
PROTO_LEN = 20
SEED_LEN = 12          # PAM-proximal nucleotides that dominate binding
RE_VALID = re.compile(r"^[ACGT]+$")


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True,
                    help="Target sequences, headers as {genome}|...")
    ap.add_argument("--metadata",
                    default="../02_genomes_raw/Clinical_Enterobacterales_discovery_8000.tsv")
    ap.add_argument("--out", default="../11_guides/candidates.tsv")
    ap.add_argument("--min-coverage", type=float, default=0.90,
                    help="Minimum cov_exact to report")
    ap.add_argument("--gc-min", type=float, default=0.35)
    ap.add_argument("--gc-max", type=float, default=0.75)
    ap.add_argument("--refs-rank", choices=["family", "genus"], default="genus",
                    help="Taxonomic rank used to stratify reference selection "
                         "and coverage reporting. GTDB collapses all of "
                         "Enterobacterales into f__Enterobacteriaceae, so "
                         "'family' yields only a handful of references and "
                         "silently undersamples the candidate space.")
    ap.add_argument("--refs-per-family", type=int, default=3,
                    help="References per group at the chosen rank")
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


def enumerate_protospacers(seq):
    """All 20mers followed by NGG, on both strands, returned 5'->3' on the
    PAM-bearing strand."""
    out = set()
    for strand_seq in (seq, revcomp(seq)):
        n = len(strand_seq)
        for i in range(n - PROTO_LEN - 3 + 1):
            proto = strand_seq[i:i + PROTO_LEN]
            pam = strand_seq[i + PROTO_LEN:i + PROTO_LEN + 3]
            if pam[1:3] == "GG" and RE_VALID.match(proto) and pam[0] in "ACGT":
                out.add(proto)
    return out


def has_match_with_pam(fwd, rev, motif):
    """True if motif occurs on either strand followed by NGG."""
    L = len(motif)
    for s in (fwd, rev):
        start = 0
        while True:
            i = s.find(motif, start)
            if i == -1:
                break
            j = i + L
            if j + 3 <= len(s) and s[j + 1] == "G" and s[j + 2] == "G":
                return True
            start = i + 1
    return False


def gc(s):
    return (s.count("G") + s.count("C")) / len(s) if s else 0.0


def main():
    args = parse_args()
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Load target sequences, one entry per genome
    # ------------------------------------------------------------------
    seqs = {}
    for header, seq in read_fasta(args.fasta):
        genome = header.split("|")[0].split()[0]
        if genome not in seqs:
            seqs[genome] = seq
    print(f"Target sequences: {len(seqs)} genomes")

    # ------------------------------------------------------------------
    # Taxonomy
    # ------------------------------------------------------------------
    fam_of, gen_of = {}, {}
    meta_path = Path(args.metadata)
    if meta_path.exists():
        meta = pd.read_csv(meta_path, sep="\t", low_memory=False)
        for _, r in meta.iterrows():
            acc = r.get("ncbi_genbank_assembly_accession")
            tax = str(r.get("gtdb_taxonomy", ""))
            if not isinstance(acc, str):
                continue
            f = re.search(r"f__([^;]*)", tax)
            g = re.search(r"g__([^;]*)", tax)
            fam_of[acc] = f.group(1) if f else "unknown"
            gen_of[acc] = g.group(1) if g else "unknown"
        matched = sum(1 for g in seqs if g in fam_of)
        print(f"  with taxonomy: {matched}")
    else:
        print(f"WARNING: {meta_path} not found; no per-family breakdown")

    rank_of = gen_of if args.refs_rank == "genus" else fam_of
    families = defaultdict(list)
    for g in seqs:
        families[rank_of.get(g, "unknown")].append(g)
    print(f"  {args.refs_rank} groups represented: {len(families)}")
    for f, members in sorted(families.items(), key=lambda x: -len(x[1])):
        print(f"    {f:30s} {len(members)}")

    # ------------------------------------------------------------------
    # Candidate enumeration from family-spanning references
    # ------------------------------------------------------------------
    refs = []
    for f, members in families.items():
        refs.extend(sorted(members)[:args.refs_per_family])
    print(f"\nEnumerating candidates from {len(refs)} reference sequences ...")

    candidates = set()
    for g in refs:
        candidates |= enumerate_protospacers(seqs[g])
    print(f"  unique candidate protospacers: {len(candidates)}")

    # Cheap filters before the expensive scan
    filtered = []
    for c in candidates:
        if "TTTT" in c:            # terminates pol III if you express the sgRNA
            continue
        if not (args.gc_min <= gc(c) <= args.gc_max):
            continue
        if re.search(r"(A{5,}|C{5,}|G{5,}|T{5,})", c):
            continue
        filtered.append(c)
    print(f"  after GC / homopolymer / polyT filters: {len(filtered)}")

    # ------------------------------------------------------------------
    # Coverage scan
    # ------------------------------------------------------------------
    print("\nScanning coverage ...")
    prepared = {g: (s, revcomp(s)) for g, s in seqs.items()}
    n_genomes = len(prepared)

    rows = []
    for k, proto in enumerate(filtered, 1):
        seed = proto[-SEED_LEN:]

        hit_exact, hit_seed = [], []
        for g, (fwd, rev) in prepared.items():
            if has_match_with_pam(fwd, rev, proto):
                hit_exact.append(g)
                hit_seed.append(g)
            elif has_match_with_pam(fwd, rev, seed):
                hit_seed.append(g)

        cov_exact = len(hit_exact) / n_genomes
        if cov_exact < args.min_coverage:
            continue

        row = {
            "Protospacer": proto,
            "Seed": seed,
            "GC": round(gc(proto), 3),
            "cov_exact": round(cov_exact, 4),
            "cov_seed": round(len(hit_seed) / n_genomes, 4),
            "n_exact": len(hit_exact),
            "n_total": n_genomes,
        }

        # Per-family exact coverage - a guide at 95% overall can still be
        # completely absent from one family, which breaks an order-level claim.
        exact_set = set(hit_exact)
        worst_fam, worst_val = None, 1.0
        for f, members in families.items():
            v = sum(1 for m in members if m in exact_set) / len(members)
            row[f"fam_{f}"] = round(v, 3)
            if v < worst_val:
                worst_fam, worst_val = f, v
        row["worst_family"] = worst_fam
        row["worst_family_cov"] = round(worst_val, 3)

        rows.append(row)

        if k % 50 == 0:
            print(f"  {k}/{len(filtered)} scanned, {len(rows)} passing")

    if not rows:
        print("\nNo candidate reached the coverage threshold.")
        print("This is a normal result at order level. Options:")
        print("  - lower --min-coverage and design a multiplexed panel")
        print("  - rerun on 23S (longer, more semi-conserved regions)")
        print("  - accept family-stratified guides instead of one universal guide")
        return

    df = pd.DataFrame(rows).sort_values(
        ["worst_family_cov", "cov_exact"], ascending=False)

    lead = ["Protospacer", "Seed", "GC", "cov_exact", "cov_seed",
            "worst_family", "worst_family_cov", "n_exact", "n_total"]
    df = df[lead + [c for c in df.columns if c not in lead]]
    df.to_csv(out_path, sep="\t", index=False)

    fa = out_path.with_suffix(".fasta")
    with open(fa, "w") as fh:
        for i, r in enumerate(df.itertuples(), 1):
            fh.write(f">gRNA_{i:03d} cov_exact={r.cov_exact} "
                     f"worst_family={r.worst_family}:{r.worst_family_cov}\n"
                     f"{r.Protospacer}\n")

    print()
    print(f"Candidates passing: {len(df)}")
    print()
    print(df[lead].head(20).to_string(index=False))
    print()
    print("Saved:")
    print(" ", out_path)
    print(" ", fa)
    print()
    print("These are candidates on SENSITIVITY only. None is usable until it")
    print("survives 14_offtarget_screen.sh against the background set.")


if __name__ == "__main__":
    main()
