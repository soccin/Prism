# What oncoanalyser 3.0.0 actually does

The tool READMEs do not always describe what the tools do in this
pipeline. Some are out of date, some describe options that oncoanalyser
does not set, and some are wrong. This file records what was checked
against the code that ran and the files it wrote, step by step. Start
here before describing methods in a paper or report.

Tool names are explained in `OA_AlgoKey.md`.

## Status labels

Each statement is marked:

| Label | Meaning |
|---|---|
| CONFIRMED | Checked in the tool's source code at the version run, the exact command run, its log or its output files |
| README | Taken from the tool's README only; not yet checked |
| VARIANCE | The tool does something other than its README says. The README is wrong or incomplete for this run |
| OPEN | Not known yet |

Counts below are ranges over the run they were checked on (14 groups:
9 with a matched normal, 5 tumor-only), with no sample names. They show
the size of each effect; they will differ in another cohort.

## How to re-confirm anything here

Do not trust this file, the READMEs or memory. Check the run itself:

1. **Exact command run.** oncoanalyser copies every task's files into
   `out/logs/<group>/`:
   - `<task>.<group>.command.sh`: the command line with every argument.
   - `<task>.<group>.command.log`: the tool's own log, which often
     includes record counts.

   For example `out/logs/<group>/pave_somatic.<group>.command.sh`. An
   option that is not in `command.sh` was not set; the tool used its
   default.
2. **Tool versions.** `out/pipeline_info/nf_core_oncoanalyser_software_mqc_versions.yml`.
3. **Defaults.** Read the tool's source code at its release tag, not the
   README. The HMF tools are in https://github.com/hartwigmedical/hmftools,
   tag `<tool>-v<version>` (for example `pave-v1.9`). Option defaults
   are set in `<Tool>Config.java`.
4. **How oncoanalyser calls the tool.** Read
   https://github.com/nf-core/oncoanalyser at tag `3.0.0`:
   - `modules/local/<tool>/main.nf`: the command template.
   - `conf/modules.config`: extra arguments (`ext.args`). In 3.0.0 it
     sets only the log level for these tools.
5. **What a step keeps.** Count records before and after it:

   ```
   sbatch ./Prism/docs/scripts/countVcfFilters.sh out RESULT_DIR
   ```

   It counts every VCF of the chains below, for all groups, and writes:
   - `filter_counts.tsv`: records by FILTER value, per file;
   - `svtype_counts.tsv`: SV files, by SVTYPE and FILTER.

   It takes about 5 minutes. For a single file:
   `zgrep -c -v "^#" FILE.vcf.gz`, or
   `bcftools query -f '%FILTER\n' FILE.vcf.gz | sort | uniq -c`
   (`module load bcftools/1.20`). The SAGE somatic VCFs are 250 to
   650 MB; count them in a job, not on the login node.

Compare each output with its input: the number of records, and the FILTER
values. After a pipeline upgrade, check every CONFIRMED statement again:
defaults change between versions.

## 1. Somatic SNVs and indels

```
SAGE   -> sage/somatic/<T>.sage.somatic.vcf.gz   candidates past SAGE's hard filters; soft filters as FILTER labels
PAVE   -> pave/somatic/<T>.pave.somatic.vcf.gz   SAGE PASS records only; adds PON and PONGnomad labels
PURPLE -> purple/<T>.purple.somatic.vcf.gz       matched: every PAVE record. tumor-only: PASS records only
```

### Matched samples (9)

| Step | Records | What happens | Status | Evidence |
|---|---|---|---|---|
| SAGE | 0.56 to 2.8 million per sample; 0.3% to 2.8% PASS | Hard filters drop low-quality candidates before writing. Soft filters (tumor quality and support, mapping quality, edge distance, jitter, and germline evidence from the normal) label the rest. TINC recovery is on (`-run_tinc`) | CONFIRMED (labels, counts, option); README (filter rules) | Counts; `sage_somatic.*.command.sh`; SAGE README "Hard Filters", "Soft Filters", "TINC" |
| PAVE | 8,022 to 27,995: exactly SAGE's PASS records, in every sample | **Drops every record that is not PASS in SAGE's file.** About 97% to 99.7% of SAGE's records go here | CONFIRMED, VARIANCE | Counts (PAVE records = SAGE PASS records in 9 of 9). Code (`pave-v1.9`): `PaveConfig.java` defaults the filter type to PASS when no `filter_type` is given; `ChromosomeTask.java` `processVariant()` returns without writing non-PASS records. `pave_somatic.*.command.sh` sets no `filter_type`. `.command.log`: "complete for N variants" against "writing M pending variants" |
| PAVE | 27% of its records labeled (over the 9 samples) | Adds FILTER `PON` (in the panel of normals: hotspot tier 6+ samples with 5+ reads, panel tier 3+ with 3+, other 3+ samples) and `PONGnomad` (gnomAD frequency above the threshold). Writes these records too. No artefact PON (passed only in targeted mode) | CONFIRMED | Counts; `pave_somatic.*.command.sh` (`-pon_file`, `-gnomad_freq_dir`, no `-write_pass_only`, no `-pon_artefact_file`); `.command.log` "loaded PON filter"; `ChromosomeTask.java` `applyFilters()` |
| PURPLE | Same as PAVE; 5,016 to 20,642 PASS | Keeps every PAVE record with its FILTER value; adds no somatic filter. Adds VAF and copy number, biallelic status, clonality, driver annotation, the REPORTED flag | CONFIRMED | Counts (equal in 9 of 9, FILTER values identical); `purple.*.command.log` "loaded N somatic variants" |

**VARIANCE, PAVE.** The PAVE README lists `read_pass_only` ("Only process
passing variants") as an option, which implies non-PASS records are kept
by default. In PAVE 1.9 the default is the opposite: only PASS records
are processed and written. `read_pass_only` is not an option in the
code; the options are `filter_type` (ALL, PASS, PANEL) and a
process-non-PASS setting.

### Tumor-only samples (5)

Without a normal, SAGE's germline filters do not run. Germline variants
are removed only by PAVE's PON and gnomAD labels.

| Step | Records | What happens | Status | Evidence |
|---|---|---|---|---|
| SAGE | 5.8 to 6.5 million per sample; 69% to 77% PASS | No germline soft filters | CONFIRMED | Counts; FILTER values contain no germline filter |
| PAVE | 4.43 to 4.64 million | Keeps SAGE PASS records, as in matched samples, except 16 to 34 records per sample | CONFIRMED; the missing records are OPEN | Counts. PAVE's code skips records on contigs other than chr1 to 22, X and Y (`ChromosomeTask.java`); SAGE ran with `-include_mt`. Not yet checked whether the missing records are chrM |
| PAVE | 98% of its records labeled `PON` and/or `PONGnomad` (over the 5 samples) | Germline variants are labeled here | CONFIRMED | Counts |
| PURPLE | 27,666 to 162,851, all PASS | **Drops every non-PASS record**, and 0 to a few PASS records (7 in total) | CONFIRMED (counts); README for the reason | PURPLE README "Tumor only mode": SNVs and indels with PON filters, and HLA SNVs and indels, are hard filtered. Not yet checked that the 7 PASS records dropped are HLA |

**Result.** The final somatic call set is the records with FILTER = PASS
in `purple.somatic.vcf.gz`. SAGE's non-PASS candidates exist only in
`sage.somatic.vcf.gz`. In tumor-only samples a germline variant stays in
the PASS set unless it is in the panel of normals or above the gnomAD
frequency threshold. The tumor-only PASS counts (27,666 to 162,851) are
higher than the matched ones (5,016 to 20,642); residual germline
variants are the likely cause (inference, not checked).

## 2. Germline SNVs and indels (matched samples only)

```
SAGE germline  -> sage/germline/<T>.sage.germline.vcf.gz
PAVE germline  -> pave/germline/<T>.pave.germline.vcf.gz
PURPLE         -> purple/<T>.purple.germline.vcf.gz
```

| Step | Records | What happens | Status | Evidence |
|---|---|---|---|---|
| SAGE germline | 939 to 1,241 per sample; 766 to 923 PASS | Germline calling at the configured regions | CONFIRMED (counts); README (rules) | Counts |
| PAVE germline | 766 to 923: exactly SAGE's PASS records | Drops non-PASS records (same default as somatic). Annotates with the germline blacklist; no PON file; gnomAD filter switched off (`-gnomad_no_filter`) | CONFIRMED | Counts (equal in 9 of 9); `pave_germline.*.command.sh`; `.command.log` |
| PURPLE | Same as PAVE; 758 to 919 PASS | Keeps every record; labels some `LOW_VAF` and/or `LOW_TUMOR_VCN` (173 records in total over 9 samples) | CONFIRMED (counts, labels); README (rules) | Counts; PURPLE README "Genotype" |

## 3. Structural variants

```
ESVEE  -> esvee/<T>.esvee.raw.vcf.gz          all breakends, no FILTER
       -> esvee/<T>.esvee.unfiltered.vcf.gz   the same breakends with filter labels
       -> esvee/<T>.esvee.somatic.vcf.gz      somatic breakends: PASS and PON only
       -> esvee/<T>.esvee.germline.vcf.gz     germline breakends: PASS only (matched samples)
PURPLE -> purple/<T>.purple.sv.vcf.gz         see below
       -> purple/<T>.purple.sv.germline.vcf.gz
```

Counts are VCF records (breakends; a two-sided SV has two records).

### Matched samples (9)

| Step | Records | What happens | Status | Evidence |
|---|---|---|---|---|
| ESVEE | 50,302 to 139,056 breakends; 1,733 to 31,746 PASS | Filters (minimum AF, support, quality, anchor length, fragment length, strand bias, PON and others) as labels in `unfiltered.vcf.gz`; then a germline/somatic split | CONFIRMED (counts, labels); README (rules) | Counts; `esvee.*.command.log` "variants passing(N) filtered(M)" and one count per filter; ESVEE README "STEP 4: Filtering" |
| ESVEE somatic file | 263 to 2,752; 92 to 2,661 PASS | Only PASS and PON somatic breakends are written; the rest are dropped | CONFIRMED | Counts (FILTER values PASS and PON only) |
| ESVEE germline file | 34,213 to 57,253, all PASS | Only PASS germline breakends are written | CONFIRMED | Counts |
| PURPLE somatic SVs | 269 to 7,427 | Every ESVEE somatic record, plus inferred breakends (FILTER `INFERRED`, 6 to 5,360 per sample), plus germline SVs recovered as somatic (below) | CONFIRMED | Counts by SVTYPE and FILTER: PURPLE records = ESVEE somatic records + INFERRED + recovered, in 9 of 9 |
| PURPLE TINC recovery | 4 of 9 samples: 28 to 223 SVs (45 to 372 breakends) | When PURPLE measures tumor-in-normal contamination, it moves germline SVs that the contamination explains into the somatic file and out of the germline file | CONFIRMED, VARIANCE | `purple.*.command.log` "tincLevel(X) recovered N germline SVs" in exactly those 4 samples; the somatic file gains and the germline file loses the same number of records |
| PURPLE germline SVs | ESVEE germline minus the recovered SVs | | CONFIRMED | Counts |

**VARIANCE, PURPLE inferred breakends.** The PURPLE README calls them
type `INF`. In the VCF they have SVTYPE `BND` and FILTER `INFERRED`, so
they are not PASS. Select them by FILTER, not by SVTYPE.

**VARIANCE, PURPLE TINC recovery of SVs.** The PURPLE README says the SV
VCF "will contain all (filtered and unfiltered) entries from the input
structural variant VCF". It does not say that PURPLE also adds germline
SVs to it, and removes them from the germline file, when it detects
tumor-in-normal contamination.

### Tumor-only samples (5)

| Step | Records | What happens | Status | Evidence |
|---|---|---|---|---|
| ESVEE | 42,213 to 52,048 breakends | Same filters; no germline file | CONFIRMED | Counts |
| ESVEE somatic file | 31,025 to 34,037; 1,390 to 3,422 PASS | PASS and PON only; most are PON | CONFIRMED | Counts |
| PURPLE somatic SVs | 1,390 to 3,422, all PASS | **Drops every PON record**, and 14 PASS records in total. **Writes no inferred breakends** | CONFIRMED (counts); README for the PON drop; OPEN for the 14 PASS records and the missing inferred breakends | Counts; PURPLE README "Tumor only mode" |

## 4. What later steps select from PURPLE's files

No later step rewrites PURPLE's VCFs. Each later tool selects records
when it reads them:

| Tool | Selection | Status |
|---|---|---|
| LINX | Leaves out PON-filtered SVs. Marks likely artefacts (short foldback inversions, isolated translocations and single breakends without copy number support, pairs of inferred breakends) as ARTIFACT and excludes them from clustering; they stay in its tables | README |
| CHORD | PASS variants only, unless `-include_non_pass` is set | README; flag not yet checked in `chord.*.command.sh` |
| LILAC | PASS variants in coding or splice regions of the HLA genes | README |
| NEO | PASS point mutations with missense, frameshift, inframe or stop-lost effect | README |
| CUPPA | Ignores passenger indels in repeats longer than 6 bases, and in MSI samples | README |
| ORANGE | PASS variants with the REPORTED flag; reported fusions and breakends; drivers not NOT_REPORTED; viruses not blacklisted; UGT1A1 removed | CONFIRMED (source code at `orange-v5.0.1`; counts match in all 14 samples). See `dbScheme/orange.md` |

## 5. Other variances from the READMEs found so far

Details are in `dbScheme/group_a.md` and `dbScheme/orange.md`.

- PURPLE `purity.tsv` `status` holds the fit method, not the QC status.
- PURPLE README: the TINC level is "measured by Sage". No SAGE output
  file holds it; it first appears in `purple.qc`.
- TEAL output columns are named differently from its README.
- LILAC `lilac.qc.tsv` has one row per gene group, not per-gene columns.
- CUPPA `vis_data.tsv` `feat_value` holds the signature count in
  `sig_quantile` rows; the README says it is empty there.
- ORANGE README: the JSON holds every mutation. It holds only the
  REPORTED ones (under 0.1% of PASS somatic variants).
- ORANGE errors in its own JSON (mislabeled CUPPA SV count, `TABBP` gene
  symbol, PARTIAL for whole-gene deletions): `dbScheme/orange.md`.

## 6. Open

- Tumor-only: the 16 to 34 SAGE PASS records per sample missing from
  PAVE's file (chrM suspected), the 7 PASS SNVs/indels and 14 PASS SVs
  PURPLE drops, and why PURPLE writes no inferred breakends.
- PURPLE TINC recovery: the rule that selects which germline SVs are
  moved (PURPLE code at `purple-v4.4`).
- SAGE hard-filter settings actually used (`sage_*.command.sh`, SAGE
  code at `sage-v5.0.2`).
- Check the `command.sh` of LINX, CHORD, LILAC, NEO, CUPPA and SIGS for
  options that change the README defaults.
