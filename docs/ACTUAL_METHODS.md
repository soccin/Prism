# What oncoanalyser 3.0.0 actually does

The tool READMEs do not always describe what the tools do in this
pipeline. Some are out of date, some describe options that oncoanalyser
does not set, and some are wrong. This file records what was checked
against the code that ran and the files it wrote, step by step. Start
here before describing methods in a paper or report.

Tool names are explained in `OA_AlgoKey.md`. What Prism does after the
pipeline (`postOA.sh` and later steps) is in `prism_methods.md`.

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
| PAVE | 27% of its records labeled (over the 9 samples) | Adds FILTER `PON` (in the panel of normals: hotspot tier 6+ samples with 5+ reads, panel tier 3+ with 3+, other 3+ samples) and `PONGnomad` (gnomAD frequency above the threshold). Writes these records too. No artifact PON (passed only in targeted mode) | CONFIRMED | Counts; `pave_somatic.*.command.sh` (`-pon_file`, `-gnomad_freq_dir`, no `-write_pass_only`, no `-pon_artefact_file`); `.command.log` "loaded PON filter"; `ChromosomeTask.java` `applyFilters()` |
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

**ESVEE `breakend.tsv`.** Its `Id` is the VCF ID of the breakend in the
ESVEE VCFs and in PURPLE's SV VCF, with the same chromosome, position
and orientation; orientation is 1 when the VCF ALT starts with the
reference base, else -1. Inferred breakends are PURPLE's and have no
row. CONFIRMED: one group fully (2,128 records; all 1,910 non-inferred
match); `scripts/svCohort.R` checks every group on every run and stops
on a mismatch.

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
| LINX | Somatic: analyzes the SVs with FILTER PASS or INFERRED and leaves out PON. Germline: keeps PASS and PON. Marks likely artifacts as ARTIFACT; they stay in its tables (details below) | CONFIRMED (selection, tables, artifact rules; the README's artifact rules differ from the code, see below) |
| CHORD | PASS variants only, unless `-include_non_pass` is set | README; flag not yet checked in `chord.*.command.sh` |
| LILAC | PASS variants in coding or splice regions of the HLA genes | README |
| NEO | PASS point mutations with missense, frameshift, inframe or stop-lost effect | README |
| CUPPA | Ignores passenger indels in repeats longer than 6 bases, and in MSI samples | README |
| ORANGE | PASS variants with the REPORTED flag; reported fusions and breakends; drivers not NOT_REPORTED; viruses not blacklisted; UGT1A1 removed | CONFIRMED (source code at `orange-v5.0.1`; counts match in all 14 samples). See `dbScheme/orange.md` |

### LINX 2.3.1

LINX makes no SV calls. It reads PURPLE's SV VCFs and copy number and
annotates, clusters and chains the SVs. `postOA.sh` aggregates its
tables and joins them to the SV events (`prism_methods.md` section 1).
Recheck with:

```
Rscript --vanilla Prism/docs/scripts/checkLinx.R out
```

| What | Status | Evidence |
|---|---|---|
| Arguments. Somatic: `-sv_vcf purple/<T>.purple.sv.vcf.gz -purple_dir purple -ref_genome_version 38 -ensembl_data_dir -known_fusion_file -driver_gene_panel -write_vis_data -write_neo_epitopes`. Germline: `-sv_vcf <T>.purple.sv.germline.vcf.gz -germline -driver_gene_panel`. No filter options | CONFIRMED | `linx_somatic.*.command.sh` (14 groups), `linx_germline.*.command.sh` (9 matched) |
| Somatic input: LINX loads every record of `purple.sv.vcf.gz` (log "loaded N SV data records"), then analyzes only the SVs with FILTER PASS, INFERRED or empty (log "procesed M SVs") | CONFIRMED | `SvFileLoader.loadSomaticVariantsFromVcf` (`AlwaysPassFilter`) and `createSvData` at `linx-v2.3.1`; log lines |
| `linx.svs.tsv` holds exactly the PASS and INFERRED SVs: 52 to 6,389 per sample. The PON records (91 to 556 per matched sample; none are left in tumor-only PURPLE files) are not in it | CONFIRMED | 14 of 14: rows = "procesed" count; 2 x two-breakend SVs + single-breakend SVs = PASS + INFERRED records |
| Inferred breakends have type `INF` in `linx.svs.tsv` | CONFIRMED | 9 of 9 matched: INF rows = FILTER INFERRED records |
| For an SV whose two breakends have the same position, `vcfIdEnd` often repeats `vcfIdStart`, so the mate's VCF ID is not in the table: 10 to 17 SVs in each tumor-only sample (65 in total). The 2 such SVs in matched samples have both IDs. Join LINX to the VCF on `vcfIdStart` | CONFIRMED | Counts |
| Germline input: LINX keeps the records with FILTER PASS or PON (or empty). `linx.germline.svs.tsv` holds every SV of `purple.sv.germline.vcf.gz` (all PASS here): 19,412 to 34,172 per sample | CONFIRMED | `GermlineFilter` at `linx-v2.3.1`; 9 of 9: log "loaded" = "procesed" = rows |
| SVs in clusters of category ARTIFACT stay in `linx.svs.tsv`: 9 to 2,302 per sample | CONFIRMED | `linx.clusters.tsv` `category` |
| ARTIFACT rules. Category ARTIFACT = resolved type LOW_VAF, DUP_BE, PAIR_INF or SGL_MAPPED_INF, assigned to single SVs before clustering; each such SV gets a cluster of its own and is kept out of further clustering and chaining. LOW_VAF for a BND or SGL: junction copy number maximum (`jcnMax`, = `svs.tsv` `junctionCopyNumberMax`) below 0.75 (INS never), no poly-A/T LINE motif at a breakend, inserted-sequence repeat class not `LINE/L1`, and every breakend at least 5,000 bases from the neighboring breakends on both sides. LOW_VAF for an INV: its two breakends neighbors on the chromosome (no breakend between) and at most 100 bases apart, and `jcnMax` below 0.75 or VAF (copy number change / copy number) at the start breakend below 0.05. DUP_BE: an SGL at most 1 base from the next breakend with the same orientation, or a two-breakend SV with an inserted sequence of at least 30 bases whose breakends each lie at most 1 base from a breakend of another SV, those two SVs being linked by assembly at their other breakends (a duplicate of an assembled templated insertion). PAIR_INF: neighboring INF breakends with opposite orientations and matching junction copy number. SGL_MAPPED_INF: an SGL with an inserted-sequence alignment within 1,000 bases of an INF breakend (same orientation) is merged with it; both get the type (`ClusteringPrep.linkSglMappedInferreds`). `jcnMax` combines the copy number change at the breakends with a read-count estimate from the tumor fragment count (VCF `VF`) (`CnJcnCalcs.calcAdjustedJcnValues`), so LOW_VAF uses both copy number and read support. This run, PASS events: LOW_VAF 4,221, DUP_BE 12; every ARTIFACT cluster has one SV; all 4,197 LOW_VAF BND/SGL have `junctionCopyNumberMax` below 0.75 (max 0.7498) and no breakend within 5,000 bases (min gap 5,065); all 24 LOW_VAF INV are at most 94 bases | CONFIRMED; README differs (5) | `ClusterClassification.getClusterCategory`, `isFilteredResolvedType`; `SvFiltering.filterBreakends`, `isIsolatedLowVafBnd`, `isLowVafInversion`, `checkCandidateSpanningVariant`; `SimpleClustering.hasLowJcn` (`LOW_JCN_THRESHOLD` 0.75); `ClusterAnalyser` at `linx-v2.3.1`; counts |
| Gene search for a breakend: every gene of Hartwig's `ensembl_gene_data.csv` (39,500 genes, all biotypes) whose span contains the breakend or whose 5' end is up to 100,000 bases after it (strand-aware: forward strand GeneStart - 100,000 to GeneEnd, reverse strand GeneStart to GeneEnd + 100,000). No window on the 3' side, except configured downstream windows (region type DOWNSTREAM; 43 cohort rows). A breakend within 250,000 bases before a 3' gene of a known IG-enhancer pair gets a 250,000-base window for every gene near it (`getPreGeneDistance`). No distance to the gene's 5' end is written (`breakend.tsv` `nextSpliceDistance` is the distance to the next splice site); there is no nearest-gene annotation | CONFIRMED; README agrees on 100 kb | `EnsemblDataCache.findGeneRegions` (hmf-common), `BreakendGenePrep.setSvGeneData`, `findGeneAnnotationsBySv`; `FusionCommon.DEFAULT_PRE_GENE_PROMOTOR_DISTANCE` 100000 (`fusion_gene_distance`, not set in the command); `FusionConstants.MAX_UPSTREAM_DISTANCE_ENHANCER_KNOWN` at `linx-v2.3.1` |
| `breakend.tsv` rows per attached gene: the canonical transcript; for driver panel genes with disruption reporting, also the panel's `additionalReportedTranscripts`; any transcript in a fusion. No row for a transcript the breakend lies past the 3' end of. By default no rows for the end of an SGL (genes from its inserted-sequence alignments) unless in a fusion (4 PASS SGLs have end rows), and none for SGL/INF breakends merged as SGL_MAPPED_INF (removed from the list in `run`). Region type IG is a record LINX creates for the enhancer partner of an IG-enhancer fusion (transcript named `@IGH` etc.), not a real transcript. 64,474 of 64,996 cohort rows are canonical | CONFIRMED; README differs (5) | `FusionDisruptionAnalyser.buildTranscriptList`, `BreakendGenePrep.createBreakendTranscriptData`, `DisruptionFinder.matchesDisruptionTranscript`, `getDisruptionGeneTranscripts`, `FusionFinder.generateEnhancerTranscript` at `linx-v2.3.1`; counts |
| `svs.tsv` `geneStart`/`geneEnd`: of the genes attached above (after LINX removes transcripts that are not canonical, disruptive or usable in a fusion), those whose span contains the breakend, or with a remaining transcript starting 1 to 10,000 bases after the breakend in the direction of transcription (skipped when the breakend faces into the gene and another gene's splice acceptor lies between); joined with `;`. For an SGL, `geneEnd` comes from the alignments of its inserted sequence (VCF INFO `INSALN`): 332 PASS SGLs have one. Check: of PASS breakends with only UPSTREAM rows, 2,530 of 2,667 with a gene are within 10 kb of a listed gene's 5' end, 5,092 of 5,155 without are farther (gene-level distance) | CONFIRMED; README differs (5) | `SvVarData.getGeneInBreakend`, `PRE_TRANSCRIPT_DISTANCE` 10000; `BreakendGeneData.breakendWithinGene`; `FusionDisruptionAnalyser.annotateTranscripts`; `SglMapping`, `StructuralVariantFactory` (`INSALN`) at `linx-v2.3.1`; counts |
| `linx.driver.catalog.tsv` = the REPORTED rows of `purple.driver.catalog.somatic.tsv`, unchanged, plus LINX rows with driver DISRUPTION (all REPORTED here; 1 to 25 per sample). PURPLE's NOT_REPORTED rows are not in it | CONFIRMED | 14 of 14: full-row match |
| `disruptive` (`breakend.tsv`) is set for every gene and transcript; the driver gene panel is not used. It is true when the breakend lies in an exon or intron of the transcript, then set false for: inferred breakends (type INF); intronic breakends of LINE clusters; DEL, DUP or INS with both ends in the same intron of the transcript; a DUP alone in its cluster unless both ends are in the same transcript (with exceptions for a duplicated first exon); and, in chains, breakends that return to the same intron within 5,000 bases without crossing a splice acceptor, assembled links inside one intron when both chain ends are outside genes, and linked intronic pairs whose other ends are both outside genes. It uses position and SV structure only: no read support, quality or ARTIFACT test. 48% of the canonical disruptive breakends in this run are in ARTIFACT clusters | CONFIRMED | `BreakendTransData` (initial value), `DisruptionFinder.markTranscriptsDisruptive` and chain checks, `SvVarData.isSimpleType` at `linx-v2.3.1`; counts |
| `reportedStatus` REPORTED on a breakend: disruptive, on a gene of the driver gene panel with disruption reporting on, and on its canonical or listed transcript. This is the only place the panel enters | CONFIRMED | `DisruptionFinder.findReportableDisruptions`, `getDisruptionGeneTranscripts` at `linx-v2.3.1` |
| Keys: `breakend.svId` and `links.lowerSvId`/`upperSvId` are in `svs.svId`; `svs.clusterId` = `clusters.clusterId`; `clusterCount` sums to the `svs` rows; fusion breakend IDs are in `breakend.id`; `drivers.clusterId` is in `clusters` or -1 | CONFIRMED | 14 of 14 |

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
- LINX README, `svs.tsv` `geneStart`/`geneEnd`: "Gene(s) overlapping"
  the breakend. The code also lists genes with a transcript starting up
  to 10 kb after the breakend, and for an SGL fills `geneEnd` from the
  inserted sequence's alignments (LINX section above).
- LINX README, `breakend.tsv`: one record per canonical transcript plus
  non-canonical transcripts in a fusion. The code also writes the
  driver panel's additional transcripts for disruption-reported genes.
- LINX README, "Artefact filtering": copy number change below 0.5 at
  both breakends, VAF below 0.05 at both ends, more than 5,000 bases,
  and resolved type ARTIFACT. The code tests the junction copy number
  maximum below 0.75, the VAF at the start breakend only, at least
  5,000 bases; ARTIFACT is the category and the resolved type is
  LOW_VAF, DUP_BE, PAIR_INF or SGL_MAPPED_INF. DUP_BE and
  SGL_MAPPED_INF are not described.

## 6. Open

- Tumor-only: the 16 to 34 SAGE PASS records per sample missing from
  PAVE's file (chrM suspected), the 7 PASS SNVs/indels and 14 PASS SVs
  PURPLE drops, and why PURPLE writes no inferred breakends.
- PURPLE TINC recovery: the rule that selects which germline SVs are
  moved (PURPLE code at `purple-v4.4`).
- SAGE hard-filter settings actually used (`sage_*.command.sh`, SAGE
  code at `sage-v5.0.2`).
- Check the `command.sh` of CHORD, LILAC, NEO, CUPPA and SIGS for
  options that change the README defaults.
- LINX: why `vcfIdEnd` repeats `vcfIdStart` for SVs whose breakends
  share a position. Whether fusion finding is limited to any gene list.
  212 PASS BND/SGL events pass the LOW_VAF copy-number and distance
  tests but are neither LOW_VAF nor LINE; not traced (the poly-A/T
  motif test on the inserted sequence was not re-run).
