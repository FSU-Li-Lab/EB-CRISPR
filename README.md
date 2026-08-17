# enterobacterales-grna-design

by Xiangpeng Li, PhD. Department of Chemistry and Biochemistry, Florida State University

A reproducible pipeline for designing dCas9 guide RNA panels that label an
entire bacterial order, for use in cell-sorting-based enrichment from complex
microbial communities.

Built to enrich Enterobacterales from human gut microbiome samples, but the
pipeline is not specific to that clade — point it at a different GTDB taxon and
it will design a panel for that instead.

---

## What problem this solves

Designing a guide for a single species is easy. Designing a small panel that
binds *every member of an order* and *nothing else in a stool sample* is a
different problem, and standard CRISPR design tools do not address it:

- **The readout is occupancy, not cleavage.** For a dCas9 labelling assay the
  design variable is *binding sites per genome*, because signal scales with
  labelled sites per cell. A guide present once in every genome is worse than
  one present seven times in 90% of them.
- **Protein conservation is not nucleotide conservation.** Synonymous variation
  destroys the 20 nt window a guide needs. In this dataset, *gyrB*, *infB* and
  *recA* yielded **zero** order-wide conserved protospacers despite being
  recovered from 98% of genomes.
- **Specificity is set by the background, not the target.** A guide is only as
  specific as the negative set it was screened against, and the negative set
  must be weighted by what is actually abundant in the sample.

## The purity budget

The core design constraint, derived rather than assumed. With `a` = target
fraction of the community, `r` = fraction of target cells labelled above gate,
`e` = fraction of non-target cells labelled above gate, `P` = required purity:

```
purity = a·r / (a·r + (1−a)·e)  ≥  P

    ⟹   e  ≤  a·r·(1−P) / (P·(1−a))
```

At `a` = 0.01, `r` = 1.0, `P` = 0.95 this gives **e ≤ 5.3 × 10⁻⁴** — fewer than
one non-target cell in ~1,900 may be labelled.

Two consequences drive the design:

1. **Counter-screening must be abundance-weighted.** A cross-reacting taxon
   must sit below ~0.05% of the community. One off-target in *Bacteroides* at
   20% abundance exceeds the budget by ~400×; the same off-target in a rare
   genus is irrelevant. A flat genome list treats these as equivalent.
2. **Recovery and purity trade linearly.** Since `e` scales with `r`, raising
   the sorting gate costs recovery but buys purity proportionally, so the gate
   is a tunable parameter rather than a fixed risk.

Implemented in `19_select_panel.py`.

---

## Pipeline

Scripts run in order from `scripts/`. Each takes `--help`.

### Dataset construction
| Script | Purpose |
|---|---|
| `01_filter_gtdb_metadata.py` | select taxon + quality thresholds from GTDB metadata |
| `02_build_discovery_set.py` | representatives + genus-stratified sample |
| `03_rename_genomes.py` | flatten NCBI download to one FASTA per accession |
| `04_build_qc_summary.py` | merge CheckM2 + QUAST |
| `05_filter_qc.py` | apply thresholds, stage for dereplication |
| `06_dereplicate.sh` | dRep at 0.99 ANI |
| `07_check_composition.py` | dereplication audit — species/genera retained |

### Annotation
| Script | Purpose |
|---|---|
| `08_run_bakta.sh` | annotate representatives |
| `09_check_bakta.py` | verify completeness, flag truncated output |

### Target loci
| Script | Purpose |
|---|---|
| `10_extract_rrna.sh` | barrnap 16S/23S with length filtering |
| `11_extract_core_genes.py` | single-copy housekeeping loci from Bakta output |

### Guide discovery
| Script | Purpose |
|---|---|
| `12_scan_conserved_grna.py` | conserved protospacer + PAM discovery |
| `13_build_background_set.py` | tiered negative set, capped per order |
| `14_offtarget_screen.sh` | Cas-OFFinder, one genome per invocation |
| `15_summarise_offtargets.py` | per-guide verdicts from raw screen output |

### Panel assembly
| Script | Purpose |
|---|---|
| `16_check_guide_overlap.py` | detect guides sharing a locus (both strands) |
| `17_check_spacing.py` | inter-site spacing vs dCas9 footprint |
| `18_count_sites.py` | photon budget, per-species minima |
| `19_select_panel.py` | purity budget + greedy panel assembly |
| `20_make_order_sheet.py` | protospacer → sgRNA, IVT templates |

See [`docs/RUNBOOK.md`](docs/RUNBOOK.md) for the full command sequence with
checkpoints and expected runtimes.

---

## Install

```bash
conda env create -f environment.yml   # main pipeline
conda activate grna-design

conda create -n casoff -c conda-forge -c bioconda cas-offinder
```

Cas-OFFinder is separated because its OpenCL dependencies conflict with the
rest. Use `DEVICE=G` if you have a working OpenCL GPU runtime — roughly an
order of magnitude faster than CPU.

---

## Pitfalls this pipeline handles

Each of these silently corrupts results if unhandled. All were encountered in
practice.

**Reference selection must be stratified by genus, not family.** GTDB assigns
the whole of Enterobacterales to `f__Enterobacteriaceae`. Stratifying candidate
enumeration by family therefore samples ~3 genomes out of thousands, and misses
most of the candidate space. Genus stratification recovered 7× more candidates
and produced the entire final panel. (`12_scan_conserved_grna.py --refs-rank`)

**CheckM2 cannot see rRNA.** It scores single-copy protein markers, so a genome
whose rRNA operons were collapsed by the assembler still reports ~100%
complete. `Chromosome`-level assemblies are frequently scaffolded with N-gaps
positioned exactly where the rRNA repeats sit. Draft assemblies undercounted
rRNA-derived signal by ~2.3× in this dataset. Copy-number calibration requires
`Complete Genome` assemblies specifically. (`18_count_sites.py --closed-only`)

**Guides on opposite strands can occupy the same site.** dCas9 unwinds the
duplex, so an occupied locus is unavailable regardless of strand. Two such
guides look completely dissimilar on direct comparison — the relationship only
appears after reverse-complementing one. Counting both inflates the photon
budget. (`16_check_guide_overlap.py`)

**Seed-only screening is uninformative at database scale.** A 12 nt seed plus
PAM has a random expectation of ~200 hits per guide across a few thousand
genomes; every guide "fails". Full-length matching at ≤2 mismatches is the
discriminating test (~5.5 random hits per guide over 32k genomes, versus ~101
at 3 mismatches).

**The gate does not protect rRNA guides.** Background bacteria also carry
multicopy rRNA, so a cross-reacting rRNA guide clears the sorting gate as
easily as a real target. The copy-number margin only buys specificity for
single-copy guides.

**Genus matching must handle GTDB suffixes.** GTDB splits polyphyletic genera
as `Klebsiella_A`, `Escherichia_B` and so on. Exact string equality silently
skips every split lineage. (`02_build_discovery_set.py`)

**Shell globs fail silently above ~32k files.** `cat dir/*.txt` exceeds ARG_MAX
and, combined with a `||` fallback, truncates output to zero without erroring.
Use `find -exec`. (`14_offtarget_screen.sh`)

---

## Reference result

Applied to Enterobacterales:

| | |
|---|---|
| Discovery set | 7,203 high-quality genomes |
| After 0.99 ANI dereplication | 4,691 representatives — 464 species, 33 genera |
| Final panel | 7 guides: 4× 23S rRNA, 3× *secY* |
| Sites per genome | 31 median, 28 at 5th percentile |
| rRNA copy number recovered | 7.2–7.3 (matches rrnDB) |
| Background screened | 32,536 GTDB species representatives + GRCh38.p14 |
| Bacterial off-targets at ≤2 mm | **0** |
| Human off-targets at ≤2 mm | 6, consistent with random expectation for 3.2 Gb |

Panel sequences: [`results/panel_final.tsv`](results/panel_final.tsv).

Note that scope excludes obligate insect endosymbionts (*Erwinia_E*,
*Serratia symbiotica*, *Pantoea carbekii*) and plant-associated genera
(*Pectobacterium*, *Dickeya*), which do not occur in the gut. These are the
only genomes falling below the sorting gate.

---

## Limitations

- The background is GTDB **species representatives**, so strain-level variation
  within abundant gut genera is not captured. Re-screening against sample-derived
  metagenomes is the appropriate next step.
- In silico specificity does not establish binding affinity. Predicted sites
  require experimental confirmation.
- The abundance weights used for the purity budget are literature-derived and
  should be replaced with values measured from the target sample type.
- Guides are designed for SpCas9 NGG PAM. Other Cas variants need the PAM
  pattern changed in `12_scan_conserved_grna.py` and `14_offtarget_screen.sh`.

## Citation

If this is useful, please cite the repository. Method description forthcoming.

## License

see [LICENSE](LICENSE).
