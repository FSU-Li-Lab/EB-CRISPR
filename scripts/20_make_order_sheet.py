#!/usr/bin/env python3
"""
28_make_grna_order_sheet.py

Convert protospacers into sequences you can actually order or transcribe.

THE THREE THINGS PEOPLE GET WRONG
---------------------------------
1. The sgRNA spacer is the protospacer sequence AS WRITTEN, with T -> U. It is
   NOT the reverse complement. The spacer base-pairs with the strand OPPOSITE
   the one carrying the PAM, which means it matches the PAM-bearing strand -
   the one reported here. Reverse-complementing is the single most common
   error in guide ordering and produces a guide that binds nothing.

2. The PAM is NOT part of the guide. NGG must be present in the target DNA,
   immediately 3' of the protospacer. Including it in the sgRNA breaks the
   complex. The PAM is never synthesised.

3. T7 transcription needs a 5' G, and initiates far better with GG. If the
   spacer does not start with G, prepend one. A 21-nt spacer with an extra 5' G
   is well tolerated by SpCas9; substituting the first base instead costs you
   a match at the PAM-distal end, which for a BINDING assay is the end you can
   most afford to lose but should still avoid if prepending works.

SCAFFOLD CHOICE MATTERS FOR YOUR APPLICATION
--------------------------------------------
Two scaffolds are in common use. The original (Jinek/Mali) works fine for
cleavage. For dCas9 IMAGING and labeling, the optimised "F+E" scaffold from
Chen et al. 2013 (Cell 155:1479) - an A-U flip removing a partial pol III
terminator, plus an extended stem-loop - gives markedly better complex
stability and signal. Since Hy-SCALE reads occupancy rather than cutting,
complex stability is the whole ballgame, so F+E is the sensible default.

VERIFY BOTH SCAFFOLD SEQUENCES against the primary source before ordering.
They are reproduced here for convenience, not as an authority.

Usage:
    python3 28_make_grna_order_sheet.py \\
        --guides ../11_guides/panel_nonredundant.fasta \\
        --out ../11_guides/order_sheet.tsv
"""

import argparse
from pathlib import Path

# Verify against the primary literature before ordering.
SCAFFOLD_STANDARD = (
    "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTT"
    "GAAAAAGTGGCACCGAGTCGGTGC"
)
SCAFFOLD_FE = (
    "GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCG"
    "TTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
)
T7_PROMOTER = "TAATACGACTCACTATAG"


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--guides", required=True)
    ap.add_argument("--out", default="../11_guides/order_sheet.tsv")
    ap.add_argument("--scaffold", choices=["fe", "standard"], default="fe")
    return ap.parse_args()


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


def main():
    args = parse_args()
    guides = read_fasta(args.guides)
    scaffold = SCAFFOLD_FE if args.scaffold == "fe" else SCAFFOLD_STANDARD

    print(f"Guides: {len(guides)}")
    print(f"Scaffold: {args.scaffold} ({len(scaffold)} nt)\n")

    rows = []
    for name in sorted(guides):
        proto = guides[name]
        gc = round(100 * (proto.count("G") + proto.count("C")) / len(proto), 1)

        needs_g = not proto.startswith("G")
        spacer_dna = ("G" + proto) if needs_g else proto
        spacer_rna = spacer_dna.replace("T", "U")

        # T7 promoter ends in G, which serves as the +1 base. If the spacer
        # already begins with G, that G is supplied by the promoter.
        ivt_body = spacer_dna[1:] if spacer_dna.startswith("G") else spacer_dna
        ivt_template = T7_PROMOTER + ivt_body + scaffold

        flags = []
        if "TTTT" in spacer_dna:
            flags.append("polyT")
        if needs_g:
            flags.append("5'G added")
        if gc < 35 or gc > 75:
            flags.append(f"GC {gc}%")

        rows.append({
            "Guide": name,
            "Protospacer_DNA_target": proto,
            "PAM": "NGG (in target, not synthesised)",
            "Spacer_DNA": spacer_dna,
            "Spacer_RNA": spacer_rna,
            "Spacer_len": len(spacer_dna),
            "GC_pct": gc,
            "IVT_template_DNA": ivt_template,
            "Flags": ";".join(flags) or "-",
        })

    try:
        import pandas as pd
        df = pd.DataFrame(rows)
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.out, sep="\t", index=False)
        print(df[["Guide", "Protospacer_DNA_target", "Spacer_RNA",
                  "Spacer_len", "GC_pct", "Flags"]].to_string(index=False))
        print(f"\nSaved: {args.out}")
    except ImportError:
        for r in rows:
            print(r)

    print()
    print("For synthesis (IDT, Synthego, etc.): order the Spacer_RNA column as")
    print("a chemically modified sgRNA, or supply Spacer_DNA and let the vendor")
    print("append their scaffold. For in-house IVT, use IVT_template_DNA.")
    print()
    print("Order NOTHING until the full background screen finishes - guides that")
    print("cross-react with abundant gut taxa will be eliminated, and synthesis")
    print("is the expensive irreversible step.")


if __name__ == "__main__":
    main()
