# Group A: one value set per sample

Group A is the per-sample summary: values that describe the sample as a
whole rather than single events. This page lists the Group A tool files
and how much of each ORANGE holds. Field-level detail is in
`tables/group_a_fields.csv` (sheet `group_a_fields` of
`dbScheme.xlsx`): every field with its description, its ORANGE path if
any, and whether it is present in matched and tumor-only samples.

## Readouts

Paths are relative to `out/<group>/`. Fields: number of columns, keys
or metric rows. In ORANGE: how many of those fields ORANGE holds.

| Readout | File | Layout | Fields | In ORANGE | Matched only |
|---|---|---|---|---|---|
| Purity fit and characteristics | `purple/<T>.purple.purity.tsv` | 1 row | 25 | 17 | no |
| PURPLE QC | `purple/<T>.purple.qc` | key-value | 13 | 10 | no |
| AMBER QC | `amber/<T>.amber.qc` | key-value | 4 | 1 | no |
| QSEE QC status | `qsee/<T>.qsee.status.tsv.gz` | 1 row per metric | 8 columns, 17 metrics | 4 metrics | no |
| BAM summary | `bamtools/<S>/<S>.bam_metric.summary.tsv` | 1 row per BAM | 29 | 0 | no |
| Read flag counts | `bamtools/<S>/<S>.bam_metric.flag_counts.tsv` | samtools flagstat text | 1 (16 lines) | 0 | no |
| HRD (CHORD) | `chord/<T>.chord.prediction.tsv` | 1 row | 8 | 6 | yes |
| Signatures | `sigs/<T>.sig.allocation.tsv` | 1 row per signature | 3 | 3 | yes |
| Tissue of origin (CUPPA) | `cuppa/<T>.cuppa.pred_summ.tsv` | 1 row per classifier | 11 | 7 | yes |
| Telomere length (TEAL) | `teal/<S>.teal.tellength.tsv` | 1 row per BAM | 13 | 0 | no |
| HLA alleles (LILAC) | `lilac/<T>.lilac.tsv` | 1 row per allele | 20 | 11 | no |
| HLA QC (LILAC) | `lilac/<T>.lilac.qc.tsv` | 1 row per gene group | 25 | 1 | no |
| Pharmacogenomic haplotypes (PEACH) | `peach/<N>.peach.haplotypes.best.tsv` | 1 row per haplotype | 6 | 6, DPYD only | yes |
| Pharmacogenomic QC (PEACH) | `peach/<N>.peach.qc.tsv` | 1 row per gene | 2 | 0 | yes |
| Viruses | `virusinterpreter/<T>.virus.annotated.tsv` | 1 row per virus | 11 | 9 | no |
| Immune receptor loci (CIDER) | `cider/<T>.cider.locus_stats.tsv` | 1 row per locus | 6 | 0 | no |
| Total | | | 202 | 75 | |

## Is all of Group A in ORANGE?

No. ORANGE holds 75 of the 202 fields.

Not in ORANGE at all:

- BAM metrics and read flag counts (coverage, duplicates, mapping and
  base quality, depth fractions)
- Telomere length (TEAL)
- Immune receptor loci (CIDER)
- The QSEE status table (ORANGE has only the path of the QSEE plot)
- PEACH QC

Partly in ORANGE:

- PURPLE purity file: no `gender`, `score`, `normFactor`,
  `diploidProportion` (and its range), `polyclonalProportion`,
  `somaticPenalty`. `runMode` is only implied (null `referenceId`).
- PURPLE QC: no `AmberGender`, `CobaltGender`, `TincLevel`. ORANGE has
  no sex call anywhere.
- AMBER QC: only the contamination rate; no QC status, consanguinity or
  uniparental disomy.
- CHORD: no remarks.
- CUPPA: no `extra_info` (top signatures and top features).
- LILAC: total fragment counts only (no unique, shared, wildcard
  counts); of the QC file only the status, copied to each allele.
- PEACH: DPYD only; UGT1A1 is dropped.
- Viruses: no `taxid`, no `blacklisted`.

In ORANGE but in no Group A file:

- Immune escape flags (computed by ORANGE)
- Signature etiology
- The CUPPA probability of every cancer type for every classifier, and
  the best call (the tool's summary file keeps the top 3 per
  classifier; the full list is in `cuppa.vis_data.tsv`)
- Four CUPPA SV feature counts (from `cuppa_data.tsv.gz`)

## Tumor-only samples

CHORD, signatures, CUPPA and PEACH need a matched normal; their files
do not exist for tumor-only samples. In the other files every field is
present, but some values cannot be computed without a normal (AMBER
contamination, the normal's BAM metrics and telomere length).

## Differences from the tool READMEs

Found while checking the curated table against real files:

- `purple.purity.tsv` `status` holds the fit method (NORMAL, SOMATIC
  and so on, same as `Method` in `purple.qc`), not the QC status the
  README describes.
- `purple.qc` `LohPercent` and `TincLevel` are not in the README output
  table. `LohPercent` holds a fraction.
- TEAL column names differ from the README: `rawTelomereLength`
  (README `sampleMixLength`), `finalTelomereLength` (README
  `telomereLength`); `sampleId`, `type` and `totalTelomericReads` are
  not in the README.
- `lilac.qc.tsv` has one row per gene group, with per-gene values as
  `gene=value` pairs in `LowCoverageBases` and `geneTypes`; the README
  describes per-gene columns. `lilac.tsv` has a `Genes` column the
  README does not list.
- bam-tools and SIGS describe none of their output columns; the
  descriptions here are from the column names and from matching values
  (QSEE `COVERAGE_ABOVE_10` equals `DepthCoverage_10`).
- Columns named `Percent` in `bam_metric.summary.tsv` and `percent` in
  `sig.allocation.tsv` hold fractions.
