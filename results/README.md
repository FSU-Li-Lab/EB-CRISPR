# results

Final panel from the Enterobacterales run described in the top-level README.

- `panel_final.tsv` — the seven guides with locus, GC and mean sites per genome
- `order_sheet_FINAL.tsv` — synthesis-ready sequences

For synthetic sgRNA, order the `Protospacer_DNA_target` column (20 nt); the
vendor appends the scaffold. The 5'G in `Spacer_DNA` is required only for T7
in vitro transcription, in which case use `IVT_template_DNA`.

The PAM is never synthesised — NGG must be present in the target genome
immediately 3' of the protospacer.
