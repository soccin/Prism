# Group A: which file holds each value

Group A is the per-sample summary: values that describe the sample as a
whole rather than single events. For each tool this page gives the
file that first holds its values, the fields no other file has, and the
files that only copy them. Field-level detail is in `tables/fields.csv`
(sheet `fields` of `dbScheme.xlsx`). Roles are defined in `README.md`.

Paths are relative to `out/<group>/`.

## Per tool

### PURPLE (purity, ploidy, sex, MSI, TMB, QC)

- Primary: `purple/<T>.purple.purity.tsv` (25 fields, all primary) and
  `purple/<T>.purple.qc` (10 of 13 keys primary).
- Copies in `purple.qc`: `Purity` and `Method` (= `purity` and `status`
  of `purity.tsv`), `Contamination` (= AMBER).
- Sex: `gender` in `purity.tsv` is the sex PURPLE used. `AmberGender`
  and `CobaltGender` in `purple.qc` are the two separate estimates.
  `gender` equals `CobaltGender` in all 14 samples; `AmberGender`
  differs in 1.
- `purple.qc` is the first file with `AmberMeanDepth`, `LohPercent` and
  `TincLevel`. The PURPLE README says SAGE measures the TINC level, but
  no SAGE output file holds it.
- The fit columns of `purity.tsv` equal row 1 of
  `purple.purity.range.tsv`, which lists every candidate fit.

### AMBER (contamination, consanguinity)

- Primary: `amber/<T>.amber.qc`, 4 keys.
- `Contamination` is copied by `purple.qc`, QSEE and ORANGE.

### bam-tools (read counts, coverage)

- Primary: `bamtools/<S>/<S>.bam_metric.summary.tsv` (29 fields) and
  `bam_metric.flag_counts.tsv` (16 samtools flagstat lines).
- Copies in the flag counts: `in total` and `duplicates` (=
  `TotalReads` and `DuplicateReads` of the summary). The QC-failed read
  counts are only in the flag counts file.
- The six `...Percent` and `CappedCoverage` columns are fractions of
  bases excluded from the coverage count (Picard `PCT_EXC_*`), not read
  fractions: `DuplicatePercent` is not `DuplicateReads / TotalReads`.
- `DepthCoverage_N` are the cumulative fractions of
  `bam_metric.coverage.tsv`.

### QSEE (QC verdict): secondary

- `qsee/<T>.qsee.status.tsv.gz` holds no primary value. Its 17 metric
  values:
  - 13 are copies of bam-tools, AMBER and PURPLE values (rounded to 4
    decimals).
  - 4 are derived: `MAPPED_PROPORTION` (flag counts mapped / total),
    `DUPLICATE_READS` (DuplicateReads / TotalReads),
    `DUAL_STRAND_READS` (DualStrandReads / TotalReads) and
    `COVERAGE_ABOVE_250` (from `bam_metric.coverage.tsv`).
- Only QSEE has `QcStatus` (PASS, WARN, FAIL) and `FailCondition` (its
  threshold).
- `qsee.vis.data.tsv.gz` repackages bam-tools, COBALT, REDUX and ESVEE
  distributions.

### CHORD (HRD), matched only

- Primary: `chord/<T>.chord.prediction.tsv`, 7 fields (plus sample ID).

### SIGS (mutational signatures), matched only

- Primary: `sigs/<T>.sig.allocation.tsv`: allocation and fraction per
  signature.

### CUPPA (tissue of origin), matched only

- Primary: `cuppa/<T>.cuppa.vis_data.tsv`, the full result. It is a long
  table:
  - probability of every cancer type for each classifier;
  - EVENT feature values and their odds ratios;
  - signature counts and their cohort quantiles;
  - training set metrics (the same for every sample).
- Secondary: `cuppa/<T>.cuppa.pred_summ.tsv`. Its top 3 classes and
  probabilities per classifier, and its `extra_info`, are all taken
  from `vis_data` (checked in all matched samples).

### TEAL (telomere length)

- Primary: `teal/<S>.teal.tellength.tsv`, one file per BAM. 6 fields are
  primary: the raw and final telomere length and four telomeric read
  counts.
- Copies: `purity`, `ploidy` (PURPLE; 1 and 2 for the normal),
  `duplicateProportion` (bam-tools `DuplicatePercent`), `meanReadDepth`
  and `gc50ReadDepth` (COBALT `cobalt.gc.median.tsv`).

### LILAC (HLA)

- Primary: `lilac/<T>.lilac.tsv` (alleles: fragment counts, copy number,
  somatic variant counts) and `lilac/<T>.lilac.qc.tsv` (QC per gene
  group).

### PEACH (pharmacogenomics), matched only

- Primary: `peach/<N>.peach.haplotypes.best.tsv` (haplotype and count
  for DPYD and UGT1A1) and `peach/<N>.peach.qc.tsv`.
- `function`, `linkedDrugs` and `prescriptionUrls` are annotations from
  PEACH's haplotype table.
- `peach.haplotypes.all.tsv` (candidate combinations) and
  `peach.events.tsv` (variants used) are not described here.

### VIRUSBreakend and Virus Interpreter (viruses)

- Primary: `virusbreakend/<T>.virusbreakend.vcf.summary.tsv`, one row
  per virus. It holds the detection data: taxonomy, kraken2 read counts,
  coverage, depth, base and mapping quality, integrations, QC status.
- Secondary: `virusinterpreter/<T>.virus.annotated.tsv`, one row per row
  of the summary.
  - Copies: `taxid`, `qcStatus` (empty becomes NO_ABNORMALITIES),
    `integrations`, `percentageCovered`, `meanCoverage`.
  - Annotations: `name` (its own taxonomy table; differs from
    VIRUSBreakend's names in most rows), `interpretation`,
    `driverLikelihood`, `blacklisted`.
  - Derived: `expectedClonalCoverage` (= tumor `MeanCoverage` x purity
    / ploidy, checked in all 13 samples with a virus) and `reported`.
  - Only this file has the reporting decision.

### CIDER (immune receptor loci)

- Primary: `cider/<T>.cider.locus_stats.tsv`, one row per IG or TCR
  locus.

### ORANGE: secondary

See `orange.md`.

## Counts

Over the 18 files described field by field (252 fields):

| Role | Fields |
|---|---|
| primary | 173 |
| copy | 29 |
| derived | 14 |
| annotation | 12 |
| key | 24 |

21 of the 29 copies are checked by the build in every sample; the other
8 (Virus Interpreter, COBALT sources) were checked once by separate
scripts.

## Tumor-only samples

CHORD, SIGS, CUPPA and PEACH need a matched normal; their files do not
exist for tumor-only samples. TEAL and bam-tools write only the tumor
file. In the other files every field is present, but some values cannot
be computed without a normal (AMBER contamination).

## Differences from the tool READMEs

- `purple.purity.tsv` `status` holds the fit method (NORMAL, SOMATIC
  and so on), not the QC status the README describes.
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
- `cuppa.vis_data.tsv` `feat_value` holds the signature SNV count in
  `sig_quantile` rows; the README says it is empty there.
- bam-tools and SIGS describe none of their output columns. Columns
  named `Percent` in `bam_metric.summary.tsv` and `percent` in
  `sig.allocation.tsv` hold fractions.
