# SV report: `scripts/svReport.R`

The SV report is the workbook people read about the somatic structural
variants (SVs) of a cohort: the events, the genes LINX puts at their
breakends, the fusions of two genes LINX calls, and for each event five flags of
evidence that it may not be a real tumor event. Events with no flag are
the high-confidence set. Gene and fusion sheets count tumors over the
high-confidence events and over all PASS events side by side.

Every statement about a Hartwig tool below was checked in its source
code at the tag the run used (LINX `linx-v2.3.1`, hmftools repository)
and, where possible, against this run's files. Where the tool's README
says something else, this is marked **README differs**. Interpretations
that are the user's, not the tool's, are marked as user decisions.
This document was checked against the code by an independent review on
2026-11-23 and corrected where it differed. Numbers are from the run of
2026-11-23 (SV table of `postOA.sh` job 13376226).

Contents:

1. Why the report exists
2. Running it
3. The workbook
4. Which events
5. Genes and fusions: what LINX provides
6. High-confidence events: the five flags
7. The genes and fusions_real sheets
8. How the SV table (Parquet) is built
9. Checks
10. Issues we had making it
11. Results of this run
12. Open points

## 1. Why the report exists

The request for the cohort analysis asked for two things about SVs:

1. summaries of how often each gene is hit, by fusion, disruption or
   other event, sortable by frequency; and
2. the genomic addresses of the breakpoints behind recurrent
   disruptions across samples, to decide whether they are artifacts.
   The request named genes that came up as frequent targets (ERBB4,
   CRLF2, CUX1).

LINX answers the first part only in part, and the second not at all:
it annotates each SV, but nothing in its output says whether the same
break is in another patient, in the patient's own normal, in a
population SV call set, or in a region where reads map poorly. Earlier
in this project we also found patient-specific artifact translocations
that the SV panel of normals (applied by ESVEE) and the callers'
thresholds do not remove.

So there are two files, with two jobs:

| File | Job | Who reads it |
|---|---|---|
| `res/prism/db/sv_cohort_filtered.parquet` (the SV table, 458 columns, `postOA.sh`) | Evidence store: every PURPLE SV record with everything PURPLE, ESVEE and LINX wrote for it, plus the evidence Prism computes (other patients, BAM recount, population sets, region tracks) | Analysts; this script |
| `res/prism/sv/prism_sv_report.xlsx` (this report) | The answer: events, genes, fusions, with a stated confidence rule | Biologists, clinicians |

The report reads only the SV table. Its gene and fusion columns could
have been taken from the LINX files directly; its flags could not,
because four of the five flags use evidence that only Prism computes.

## 2. Running it

```
Rscript --vanilla Prism/scripts/svReport.R [DB_DIR] [OUT_DIR]
```

Run from the project folder, after `postOA.sh` (section 8). `DB_DIR`
defaults to `res/prism/db`, `OUT_DIR` to `res/prism/sv`. It takes about
a minute and little memory. A rerun replaces both outputs. Rerun it
after every rebuild of the SV table.

Outputs:

- `res/prism/sv/prism_sv_report.xlsx`: the workbook (section 3).
- `res/prism/sv/prism_sv_events.csv`: the `events_all` sheet as CSV.

It needs `Rscript` with the R packages `tidyverse`, `fs`, `nanoparquet`
and `openxlsx`. The workbook is written with an explicit window size,
Arial font, a frozen header row on every sheet and filters on the four
data sheets.

## 3. The workbook

| Sheet | One row is | Rows (this run) |
|---|---|---|
| `fusions_real` | a LINX fusion name over the cohort, fusions of two different genes only | 1,066 |
| `genes` | a gene named at a breakend, disrupted, or in a fusion, in 2 or more tumors by at least one count | 2,435 |
| `events_hc` | a high-confidence event (no flag) | 4,180 |
| `events_all` | a PASS somatic SV event, with its flags | 12,216 |
| `rule` | a flag, term or exclusion, with its definition and count | 9 |
| `glossary` | a column name, the sheets that have it, and what it holds | 57 |

Three choices made by the user on 2026-11-23 shape these sheets:

- **Same-gene fusions are not counted as fusions.** A LINX fusion whose
  5' and 3' breakends are in the same gene (`gene_5p` = `gene_3p`) joins
  two parts of one gene, not two genes; the user treats these as
  cryptic splice products, not fusions. (LINX itself calls no splicing:
  it has already dropped same-gene pairs of different transcripts,
  non-coding pairs and pairs within one intron,
  `FusionFinder.isIrrelevantSameGene`; what remains joins two exons or
  introns of one transcript through an SV, an intragenic
  rearrangement. This run: 96 fusions from 70 DEL, 27 DUP and 4 INV
  events, none reported, all `reportedType` NONE.) They are left out of
  `fusions_real` and of the fusion counts of `genes`, and kept in the
  `fusions` and `reported_fusions` columns of the events sheets, which
  list every LINX fusion of the event. This run: 96 such fusions on 101
  PASS events in 12 tumors (row `same_gene_fusions` of `rule`). These
  genes still count as disrupted and at a breakend.
- **`genes` lists only genes with 2 or more tumors** in at least one of
  its six counts.
- **`fusions_real` and `genes` have a `samples` column** (before
  `hc_events`): the tumors behind the all-PASS counts.

In the events sheets the six gene columns are 12 characters wide; the
other columns are sized to their content.

The glossary describes every column of every sheet once (a column on
several sheets has one row listing them, joined with ", " in sheet
order). The script stops if a column has no description or a
description matches no column, so no column can be added or removed
without the glossary; whether a description's text is right is checked
by reading, as in this document.

In the other sheets, cells with several values (genes, fusions, event
keys, groups, samples) hold them sorted and joined with `;`. An empty
cell means none, except where a test does not apply: for example
`other_patients_with_event` is empty for every one-breakend event, and
`end_*` columns are empty for single breakends unless LINX filled them
(section 5.3).

## 4. Which events

One row per SV event of PURPLE's somatic SV VCF with FILTER `PASS`. In
the SV table the FILTER values are:

| FILTER | Events | In the report |
|---|---|---|
| PASS | 12,216 | yes |
| INFERRED | 8,241 | no: copy-number breakpoints PURPLE inferred from read depth, not SV calls |
| PON | 1,433 | no: matched the SV panel of normals, a filter ESVEE applies (`-pon_sv_file`, `-pon_sgl_file` in `esvee.*.command.sh`); PURPLE keeps the label; LINX does not analyze them |

The table has no event with more than one FILTER value. LINX analyzed
every PASS event (the script stops otherwise). An event is one row of
the SV table: a two-breakend SV (DEL, DUP, INV, INS, BND) or a single
breakend (SGL, only one side mapped). `sv_type` is LINX's type.

| sv_type | PASS events | High confidence |
|---|---|---|
| BND (translocation) | 6,249 | 1,883 |
| DEL | 2,001 | 897 |
| DUP | 1,227 | 449 |
| INS | 117 | 21 |
| INV | 551 | 382 |
| SGL | 2,071 | 548 |

`event_key` (`sample_id:VCF ID of the start breakend`) is unique over
the cohort and is the row key. LINX's own `svId` is numbered within a
sample; join to LINX files on `sample_id` and `svId` (in the SV table:
`linx_svId`).

## 5. Genes and fusions: what LINX provides

All gene and fusion columns of the report come from LINX. On the events
sheets they are unchanged apart from sorting and de-duplication; the
`fusions_real` sheet and the gene fusion counts leave out same-gene
fusions, and the `genes` sheet leaves out genes with fewer than 2
tumors in every count (section 3). No other gene annotation is used
(decision of 2026-11-23: no annotator is written here). It is
therefore essential to know exactly what LINX annotates. LINX was run
with `-ensembl_data_dir`, `-known_fusion_file` and `-driver_gene_panel`
and no option that changes the distances below
(`linx_somatic.*.command.sh`).

### 5.1 Which genes LINX attaches to a breakend

Source: `EnsemblDataCache.findGeneRegions` (hmf-common) and
`BreakendGenePrep.setSvGeneData` / `findGeneAnnotationsBySv` (LINX).

- Genes come from Hartwig's Ensembl file
  (`common/ensembl_data/ensembl_gene_data.csv`, 39,500 genes, every
  biotype: protein coding, lncRNA, pseudogenes, snRNA and others).
- A gene is attached when the breakend lies in the gene's span, or up
  to 100,000 bases before its 5' end, taking the strand into account
  (forward strand: from GeneStart - 100,000 to GeneEnd; reverse strand:
  from GeneStart to GeneEnd + 100,000). The 100,000 is LINX's
  `fusion_gene_distance` default (`FusionCommon.DEFAULT_PRE_GENE_PROMOTOR_DISTANCE`),
  used because fusion finding is on.
- There is no window on the 3' side in general: a breakend 1 base past
  a gene's 3' end is not attached to that gene. Two exceptions are
  configured in LINX:
  - Downstream windows for a few genes from the known fusion file (the
    README names BCL2 and MYC; `SpecialFusions`, then
    `EnsemblDataCache.findGeneRegions`). Breakends found that way have
    region type `DOWNSTREAM` (43 transcript rows in the cohort, 31 on
    PASS events).
  - A 250,000-base window (`MAX_UPSTREAM_DISTANCE_ENHANCER_KNOWN`) when
    the breakend lies within 250,000 bases before the 5' end of a 3'
    gene of a known IG-enhancer pair. The larger window is then used for
    every gene near that breakend, not only for that gene
    (`BreakendGenePrep.getPreGeneDistance`).
- No distance from the breakend to the gene's 5' end is written. The
  only distance in `breakend.tsv` is `nextSpliceDistance` (in the SV
  table `{side}_linx_bnd_nextSpliceDistance`): the distance to the next
  splice acceptor or donor in the direction the breakend faces. For a
  breakend before a gene and facing it, that is the distance to the
  start of the transcript's second exon (its first splice acceptor),
  not to the gene.

The README agrees on 100 kb ("within 100kb upstream of the 1st base of
the transcript").

### 5.2 `breakend.tsv`: one row per breakend and transcript

Source: `FusionDisruptionAnalyser.buildTranscriptList`,
`BreakendGenePrep.createBreakendTranscriptData`,
`DisruptionFinder.matchesDisruptionTranscript`.

For each attached gene LINX writes:

- the canonical transcript;
- for a driver gene panel gene with disruption reporting on, also the
  transcripts the panel lists for it (`additionalReportedTranscripts`;
  `DisruptionFinder.getDisruptionGeneTranscripts`);
- any other transcript used in a fusion.

A transcript gets no row when the breakend lies past its 3' end.
Single breakends: genes attached to the end of a single breakend
through the alignments of its inserted sequence (5.3) get rows only
when they are in a fusion ("by default don't include SGL alt-mappings"
in `buildTranscriptList`; 4 PASS single breakends have end rows), and
single or inferred breakends that LINX linked to each other
(SGL_MAPPED_INF, 6.4) get no rows. In the cohort 64,474 of the 64,996
rows are canonical. `regionType` is
EXONIC or INTRONIC (in the transcript), UPSTREAM (before it, within the
window of 5.1), DOWNSTREAM (special windows only), or IG: a record LINX
creates for the immunoglobulin enhancer partner of an enhancer fusion,
with a transcript named after the region (`@IGH` and the like;
`FusionFinder.generateEnhancerTranscript`), not a real transcript.
`disruptive` and `reportedStatus` are defined in `ACTUAL_METHODS.md`
(LINX section; both checked in the code).

**README differs:** it says one record per canonical transcript plus
records for non-canonical transcripts in a fusion; the code also writes
the panel-listed transcripts of driver genes.

### 5.3 `svs.tsv` `geneStart` / `geneEnd`: one gene list per breakend

Source: `SvVarData.getGeneInBreakend`, `BreakendGeneData.breakendWithinGene`
(`PRE_TRANSCRIPT_DISTANCE = 10000`).

From the genes attached in 5.1, after LINX has removed transcripts that
are neither canonical, disruptive nor usable in a fusion, a gene is
listed when:

- the breakend lies in the gene's span (GeneStart to GeneEnd, ends
  included); or
- one of its remaining transcripts starts 1 to 10,000 bases after the
  breakend in the direction of transcription (forward strand:
  TransStart - position; reverse strand: position - TransEnd). A
  transcript is skipped for this test when the breakend faces into the
  gene and a splice acceptor of another gene lies between the breakend
  and the transcript (`hasNegativePrevSpliceAcceptorDistance`).

Several genes are joined with `;` by LINX. For a single breakend (SGL)
`geneEnd` is filled from the alignments of its inserted sequence (PURPLE
SV VCF INFO `INSALN`, parsed by `SglMapping`), not from a breakend
position: 332 PASS single breakends have a `geneEnd`.

**README differs:** it says "Gene(s) overlapping start breakend of SV".
Checked on this run: of the PASS breakends whose LINX rows are all
UPSTREAM, 2,530 of the 2,667 with a `geneStart`/`geneEnd` are within
10 kb of a listed gene's 5' end, and 5,092 of the 5,155 without one are
farther (distance measured to the gene's 5' end, not the transcript's,
so a few exceptions are expected).

### 5.4 What LINX does not provide

- **No nearest gene.** LINX names genes a break could disrupt or fuse,
  not the closest gene. Of the 22,361 breakends of PASS events with a
  position: 8,240 (36.8%) are in a protein-coding gene (canonical
  transcript INTRONIC or EXONIC); 5,974 (26.7%) are only upstream (up to
  100 kb) of a protein-coding gene; 3,533 (15.8%) have only non-coding
  genes; 4,614 (20.6%) have no gene at all.
- **No distance** to the gene (only `nextSpliceDistance`, 5.1).
- **Nothing on the 3' side** of a gene, except the few configured
  downstream windows (5.1).
- **Genes wholly inside a DEL or DUP** are in no per-breakend column:
  the breakend rows cover genes at the two ends only. The one place
  they can appear is LINX's `drivers.tsv` (SV table
  `linx_cluster_driver_gene`), which links PURPLE driver genes to the
  SV cluster that explains their copy number change; this run has 9
  such links on single-SV PASS DEL or DUP clusters (for example a
  deletion with MTAP, CDKN2A and CDKN2B as drivers, none at its
  breakends). The report does not show `linx_cluster_driver_*`.
  Gene-level copy number is in PURPLE's gene copy number files.

### 5.5 The report's gene and fusion columns

| Column | From |
|---|---|
| `start_gene`, `end_gene` | `svs.tsv` `geneStart`, `geneEnd` (5.3), each gene once, sorted. Usually one gene: of the 22,361 PASS breakends, 11,087 have one, 1,587 two or more, 9,687 none (LINX lists one gene twice at one breakend) |
| `start_disrupted_gene`, `end_disrupted_gene` | `breakend.tsv` rows with `canonical` and `disruptive` TRUE |
| `start_driver_disruption`, `end_driver_disruption` | `breakend.tsv` rows with `reportedStatus` REPORTED |
| `fusions`, `reported_fusions` | `fusion.tsv` `name` (and `reported`) of every fusion the event takes part in, each once, including fusions whose 5' and 3' gene are the same (left out of `fusions_real`; section 3) |
| `role` (genes sheet) | driver gene panel `likelihoodType` (ONCO, TSG), via the SV table's `*_gene_role`, which exists only for genes with `breakend.tsv` rows. A panel gene seen only in `start_gene`/`end_gene` (for example a single breakend's insert-alignment gene) would get an empty role; there are none in this run (checked against `DriverGenePanel.38.tsv`) |

A fusion is listed once per event. In the SV table each fusion appears
once for its 5' breakend's SV and once for its 3' breakend's SV
(`linx_fusion_side`), so a fusion made by one SV appears twice on that
SV's row (1,149 fusions) and a fusion made by two chained SVs once on
each (22 fusions); see issue 5 in section 10. Fusion names are unique
within a sample in `fusion.tsv` (1,171 names, 1,171 rows), so keeping
one entry per name and event loses nothing. LINX names an enhancer
5' partner by its region (`@IGH_PSD3`); the gene at that breakend is in
the `fusions_real` sheet's `gene_5p` (from the SV table's
`fusion_gene_5p`). Same-gene fusions are recognized by these breakend
genes (`gene_5p` = `gene_3p`), not by the name, so an enhancer fusion
is never one.

## 6. High-confidence events: the five flags

An event is **high confidence when none of the five flags is set**
(user decision, 2026-11-23). A flag is evidence that the event may not
be a real tumor event; it is not proof, and a real event can carry one.
No event is removed: `events_all` holds every PASS event with its
flags, and `events_hc` is the subset with `n_flags` = 0.

| Flag | Set when | SV table columns | Computed by |
|---|---|---|---|
| `flag_own_normal` | split reads for the break in the group's own normal BAM, either end | `{start,end}_bam_normal_split_n` > 0 | Prism BAM recount (`bamEvidence.py`) |
| `flag_other_patients` | the same event in another group's calls, or split reads for it in another group's BAM | `xs_pair_somatic_groups_w50`, `xs_pair_germline_groups_w50` (two-breakend events); `start_xs_somatic_groups_w50`, `start_xs_germline_groups_w50` (one-breakend events); `{start,end}_bam_other_{normals,tumors}_split_samples` (any event) > 0 | Prism cross-sample matching (`svCohort.R`), BAM recount |
| `flag_population` | a population SV site matches both breakpoints | `{gnomad_sv,kgp_sv,dbvar_nstd166,dbvar_nstd200}_bp100_n` > 0 | Prism population matching (`svCohort.R`, `popSvSites.sh`) |
| `flag_linx_artifact` | LINX puts the event in a cluster of category ARTIFACT | `linx_cluster_category` = ARTIFACT | LINX |
| `flag_region` | a breakend in a low-mappability or blacklisted region | `{start,end}_region_{giab_lowmappability,encode_blacklist,hmf_sv_prep_blacklist,hmf_unmap}` TRUE | Prism region tracks (`svCohort.R`, `fetchTracks.sh`) |

An empty value (for example no normal) counts as not set. The exact
definitions of the SV table columns are in `prism_methods.md`
section 1.2; the essentials follow.

### 6.1 `flag_own_normal`

`bam_normal_split_n` is a recount in the group's own normal BAM: reads
clipped by at least 5 bases at the breakend (within its confidence
interval plus 5 bases, on the breakend's side) whose supplementary
alignment (`SA` tag) ends at the mate breakend by the same rule.
Duplicate, secondary, QC-fail and unmapped reads are skipped; no
mapping-quality threshold. It is empty for tumor-only groups and for
single breakends (no mate), so the flag can only be set for
two-breakend events of the 9 groups with a matched normal
(`has_normal`). Set on 34 events.

### 6.2 `flag_other_patients`

The other groups' calls are three sets per group: **somatic**, every
record of PURPLE's `purple.sv.vcf.gz` except FILTER INFERRED (so PASS
and PON records); **germline**, every record of
`purple.sv.germline.vcf.gz` (matched groups only); and unfiltered
ESVEE calls, which this flag does not use.

- Two-breakend events: another group has two records of one SV
  (linked by `MATEID`), one within 50 bases of the event's start and
  one within 50 bases of its end, each with the same orientation, in
  its somatic or germline set.
- One-breakend events: another group has a record within 50 bases of
  the breakend, with the same orientation, in those sets.
- Any event: another group's normal or tumor BAM has at least one split
  read for the break (as in 6.1), at either end. This part is empty,
  so never set, for single breakends.

Positions only are compared: no confidence interval, homology or
inserted sequence. A real event that occurs at the same position in
several tumors is flagged as well. `other_patients_with_event` lists
the groups found by the two-breakend call match only. Set on 2,338
events.

### 6.3 `flag_population`

A site of a population SV call set matches when one of its breakpoints
is within 100 bases of the event's start and its other breakpoint
within 100 bases of the event's end (same chromosome, either way
round), and its type is allowed for the event's VCF `SVTYPE`: DEL with
DEL, CNV or BND sites; DUP with DUP, CNV or BND; INV with INV or BND;
INS with INS; BND with BND or CTX. Any FILTER counts. The sets are
gnomAD-SV v4.1, the 1000 Genomes high-coverage set (3,202 samples),
gnomAD-SV v2.1 (dbVar nstd166) and CCDG (dbVar nstd200), all short-read.
One-breakend events are never flagged (the rule needs two breakpoints).
101,132 gnomAD BND records carry no second chromosome position; their
second breakpoint is taken as their `END`, which equals `POS` in 98,515
of them (`prism_methods.md`, Population SV call sets). Set on 1,006
events.

### 6.4 `flag_linx_artifact`

Source: `ClusterClassification.getClusterCategory` and
`isFilteredResolvedType`, `SvFiltering.filterBreakends`,
`SimpleClustering.hasLowJcn` (`LOW_JCN_THRESHOLD = 0.75`). LINX gives a
cluster the category ARTIFACT when its resolved type is one of four
types. Each is assigned to single SVs before clustering; the SV is put
in a cluster of its own and kept out of all further clustering and
chaining (`ClusterAnalyser`: `applyFilters`, `clusterExcludedVariants`).
In this run every ARTIFACT cluster has one SV (7,001 events with any
FILTER, 4,233 of them PASS).

| Resolved type | Rule in the code | PASS events |
|---|---|---|
| LOW_VAF, BND or SGL | LINX's junction copy number maximum (`jcnMax`, written as `svs.tsv` `junctionCopyNumberMax`) below 0.75; no poly-A/poly-T LINE motif at a breakend and inserted-sequence repeat class not `LINE/L1`; and every breakend at least 5,000 bases from the next breakend on both sides (all breakends LINX analyzed in the sample) | 3,610 BND, 587 SGL |
| LOW_VAF, INV | the inversion's two breakends are neighbors on the chromosome (no other breakend between them) and at most 100 bases apart, and `jcnMax` below 0.75 or the allele fraction at the start breakend (copy number change / copy number) below 0.05 | 24 |
| DUP_BE | a single breakend (SGL) at most 1 base from the next breakend, with the same orientation; or a two-breakend SV with an inserted sequence of at least 30 bases (`MIN_TEMPLATED_INSERTION_LENGTH`) whose two breakends each lie at most 1 base from a breakend of another SV (same orientation), where those two other SVs are linked by assembly at their other breakends, the ends of an assembled templated insertion (`checkCandidateSpanningVariant`): a duplicate call of that insertion | 12 |
| PAIR_INF | two neighboring inferred breakends with opposite orientations and matching junction copy number | 0 (inferred breakends are not PASS) |
| SGL_MAPPED_INF | a single breakend with an alignment of its inserted sequence within 1,000 bases of an inferred breakend, same orientation, is merged with it (`ClusteringPrep.linkSglMappedInferreds`); both the single and the inferred breakend get this type | 0 |

Checked on this run: all 4,197 LOW_VAF BND and SGL events have
`junctionCopyNumberMax` below 0.75 (largest 0.7498) and no other
breakend within 5,000 bases (smallest gap 5,065); all 24 LOW_VAF
inversions are at most 94 bases long. The reverse is not exact: 499
PASS BND and SGL events also pass the copy-number and distance tests
but are not LOW_VAF. 287 of them ended in LINE clusters, which agrees
with the LINE exemption but does not prove it (the exemption is tested
per SV before clustering; the LINE type is given later); the other 212
(UNBAL_TRANS 83, COMPLEX 73, PAIR_OTHER 31, SGL 20, others 5) were not
traced (the poly-A/T test reads the inserted sequence, which was not
re-tested here).

**README differs** ("Artefact filtering"): it gives copy number change
below 0.5 at both breakends and VAF below 0.05 at both ends, and more
than 5,000 bases; the code tests the junction copy number maximum below
0.75, the VAF at the start breakend only, and at least 5,000 bases. It
says the filtered variants get resolved type ARTIFACT; ARTIFACT is the
category, the resolved type is one of the four above. It does not
mention DUP_BE or SGL_MAPPED_INF.

What `jcnMax` measures: LINX estimates the junction copy number from up
to three observations, the copy number change at each breakend and a
read-count estimate (PURPLE's junction copy number scaled by a Poisson
interval on the tumor fragment count of the start breakend, VCF `VF`),
weighted by their uncertainties; `jcnMax` is the upper end of that
estimate (`SvCNData.setSvCopyNumberData`,
`CnJcnCalcs.calcAdjustedJcnValues`). So LOW_VAF is LINX's judgment that
an isolated SV has too little copy-number change and read support
combined. Set on 4,233 events (4,221 LOW_VAF, 12 DUP_BE).

### 6.5 `flag_region`

A breakend is in a track when its position lies in one of the track's
intervals: a single base, no window, after every interval is turned
into 1-based coordinates with both ends included (BED files: start + 1;
Hartwig's `unmap_regions.38.tsv` is already 1-based). Four tracks:

| Track | Source |
|---|---|
| `giab_lowmappability` | GIAB v3.6 GRCh38 low mappability (`GRCh38_lowmappabilityall.bed.gz`) |
| `encode_blacklist` | ENCODE blacklist v2 (hg38) |
| `hmf_sv_prep_blacklist` | Hartwig `dna/sv/sv_prep_blacklist.38.bed` (no task of this run reads it) |
| `hmf_unmap` | Hartwig `common/unmap_regions.38.tsv` (the file ESVEE was given) |

Downloads, versions and md5 sums: `prism_methods.md`, Coordinates and
region tracks. Set on 3,363 events.

### 6.6 What each flag removes

| Flag | Events flagged | Flagged by this flag only |
|---|---|---|
| `flag_linx_artifact` | 4,233 | 3,012 |
| `flag_region` | 3,363 | 1,334 |
| `flag_other_patients` | 2,338 | 884 |
| `flag_population` | 1,006 | 339 |
| `flag_own_normal` | 34 | 4 |
| no flag (high confidence) | 4,180 of 12,216 | |

Number of flags per event: 0 for 4,180; 1 for 5,573; 2 for 1,995; 3 for
461; 4 for 7.

### 6.7 Limits of the rule

- Not every event can get every flag. Tumor-only groups (5 of 14) have
  no own-normal test. Single breakends cannot get `flag_own_normal` or
  `flag_population`, and their `flag_other_patients` uses the call
  match only. These events face fewer tests, so "no flag" means less
  for them: high confidence is 2,149 of 6,151 events in tumor-only
  groups and 2,031 of 6,065 in matched groups.
- The other-patient and population tests compare positions within 50 or
  100 bases; they do not compare sequence.
- Only the LINX test uses the event's own tumor evidence (copy number
  change and tumor fragment count, through `jcnMax`, and only for
  isolated BND, SGL and short INV events). The other four flags look at
  the normal, other samples, population sets and the reference. The
  report shows `tumor_fragments`, `start_af` and `end_af` but sets no
  threshold on them.

## 7. The genes and fusions_real sheets

Counts are numbers of tumors: distinct `sample_id` (one tumor per group
in this cohort). Each count is made twice, over high-confidence events
(`_hc_`) and over all PASS events (`_all_`).

`genes` (one row per gene with 2 or more tumors in at least one of the
six counts), sorted by `tumors_hc_disrupted`, then `tumors_hc_fusion`,
then `tumors_hc_at_breakend`, then gene name:

| Kind | A tumor counts when one of its events |
|---|---|
| `disrupted` | has the gene in `start_disrupted_gene` or `end_disrupted_gene` |
| `fusion` | takes part in a LINX fusion of two different genes with the gene as `gene_5p` or `gene_3p` |
| `at_breakend` | has the gene in `start_gene` or `end_gene` |

`samples` lists the tumors behind any of the all-PASS counts;
`hc_events` lists the high-confidence events behind any of the counts.

`fusions_real` (one row per fusion name over the cohort, fusions of two
different genes only), sorted by `tumors_hc`, then `tumors_all`, then
name: `gene_5p` and `gene_3p` (all values seen, normally one),
`reported` (TRUE if LINX reports it in any tumor), the two tumor
counts, `samples` (the tumors with the fusion) and `hc_events`. A
fusion made by two chained SVs is on two events; it counts as high
confidence when either event is. In this run 17 of these fusions are on
two events, and for 7 of them exactly one of the two events is high
confidence.

## 8. How the SV table (Parquet) is built

The report's only input is `res/prism/db/sv_cohort_filtered.parquet`,
written by `postOA.sh`. Full definitions are in `prism_methods.md`
section 1; this is the outline.

Inputs made first (once per site, or once per run):

| Step | Writes |
|---|---|
| `scripts/fetchTracks.sh` (once per site) | region tracks (GIAB v3.6, ENCODE blacklist v2, UCSC RepeatMasker and centromeres) and the four population SV call sets, with `tracks.tsv` (URL, date, md5) |
| `scripts/popSvSites.sh` (job, once per site) | one site table per population call set (`bcftools query`, every record kept) |
| `bamEvidence.sh` (array job, one task per BAM, per run) | read counts at every breakend of the cohort in every tumor and normal BAM (`scripts/bamEvidence.py`) |

`postOA.sh` (one job) then runs four steps into `res/prism/db/`:

1. **`scripts/postOA.R`** stacks each per-sample LINX table of every
   group into `data/linx_<mode>_<table>.tsv.gz`, lines copied as LINX
   wrote them, with `group_id` and `sample_id` in front; writes
   `samples.tsv`.
2. **`scripts/svCohort.R`** builds the SV table: one row per event of
   PURPLE's somatic `purple.sv.vcf.gz`, all FILTER values, nothing left
   out. It joins every VCF INFO and FORMAT field, the ESVEE
   `breakend.tsv` row (VCF ID = ESVEE `Id`), the LINX tables (on
   `vcfIdStart`; transcripts, fusions, links and cluster drivers as
   lists), the driver gene panel role of each gene, the 5' and 3' gene
   of each fusion, the region tracks, the same breakend or pair in
   other groups' calls, the BAM counts, and the population site
   matches. It unpacks the table for analysis (fixed pairs split,
   repeated columns dropped, fields equal on both breakends stored once
   per event), checking each step on every row. It compares the column
   dictionary with the reviewed copy `docs/sv_cohort_filtered.columns.csv`;
   if they differ it writes the new one as
   `sv_cohort_filtered.columns.new.csv` and stops, writing neither the
   table nor `sv_cohort_filtered.columns.tsv` (both are written only
   after the check passes).
3. **`scripts/tsvToParquet.py`** writes the typed Parquet copy, types
   from the dictionary; any value that does not fit stops it.
4. **`docs/scripts/checkSvCohort.R`** rebuilds every PURPLE VCF record,
   ESVEE row and LINX row from the table and compares them with the
   source files, recomputes the gene roles and fusion genes, and checks
   every Parquet cell against the TSV. Any failure fails the job.

After a change to the computed columns, `docs/scripts/spotCheckSvCohort.R`
recomputes the region, cross-sample, BAM and population columns of a
random sample by a second method (bcftools, samtools) and compares.
Last run: job 13361310, 31,772 of 31,772 values equal, on the table of
`postOA.sh` job 13347029. The table the report reads (job 13376226) was
rebuilt from scratch from the same code and inputs and passed every
check of step 4; the spot check was not rerun on it.

In R, `Prism/examples/readSvCohort.R` loads the table with helpers for
the list columns (`sv_flatten()` joins every list with `|`).

## 9. Checks

`svReport.R` stops:

- if a PASS event is not in LINX, or `event_key` is not unique;
- if a fusion has no 5' or 3' gene;
- if a gene has two driver gene panel roles;
- if the events sheet does not have one row per PASS event, or a
  `high_confidence` value is missing;
- if a column of any sheet has no glossary description, or a
  description matches no column.

Checked by hand on 2026-11-23 (scripts in the session scratch folder,
not part of the repo): the LINX rules of section 5 and 6.4 against the
LINX source and this run's data, as reported there; no missing values
in the LINX transcript lists the report reads; fusion names unique per
sample in `fusion.tsv`, and each name listed on an event either twice
(both sides, same `reported` value) or once (chained); the numbers in
sections 3, 4, 6 and 7 recomputed from the written workbook.

## 10. Issues we had making it

1. **The SV table was taken for the deliverable.** It has 458 columns
   and nothing said which ones answer which question. This report and
   this document state the purpose: the table is the evidence store,
   the report the answer.
2. **No "closest gene".** The first request was the nearest gene to
   each breakend. LINX does not provide it (5.4). Decision: use LINX
   only, no new annotation; `start_gene`/`end_gene` are LINX's
   per-breakend genes and the disrupted-gene columns show which of
   them are actually broken.
3. **`geneStart`/`geneEnd` did not match the LINX README.** The README
   says "overlapping"; the code adds genes up to 10 kb upstream (5.3).
   An earlier note in this project said not to use these columns
   because they differed from the genes the break falls inside; with
   the rule known, they are used, next to the disrupted-gene columns.
4. **Single breakends have an end gene but no end position.** LINX
   takes it from the alignments of the inserted sequence (5.3).
5. **Each fusion appeared twice on one event in the SV table.**
   `svCohort.R` attaches a fusion to the SV of its 5' breakend and to
   the SV of its 3' breakend, with `linx_fusion_side`, so that both
   SVs of a chained fusion carry it and step 4 can rebuild LINX's
   fusion table. For a one-SV fusion that lists it twice on the same
   row. LINX's `fusion.tsv` has one row per fusion. The report lists
   each fusion once per event (exact: 5.5). The SV table is unchanged.
6. **Enhancer fusion names.** LINX names an enhancer 5' partner by its
   region (`@IGH`); `fusion_gene_5p` holds the gene at the breakend.
7. **List columns.** The Parquet file stores lists as arrays. The TSV
   joins each list with the separator of its source, given in the
   dictionary's `list_sep`: `|` for 110 columns, `;` for 5 LINX
   `svs.tsv` lists (including `linx_geneStart`/`linx_geneEnd`), `,` for
   2 VCF lists. `sv_flatten()` joins every list with `|`; the report
   joins with `;`.
8. **Row counts that did not match.** LINX `svs.tsv` (20,457 rows over
   the cohort) has PASS and INFERRED events only; the SV table (21,890)
   adds the PON records; the report has PASS only (12,216).
9. **LINX IDs restart in every sample.** `svId`, `clusterId` and
   breakend `id` are unique within a sample only; join on `sample_id`
   too. `event_key` is unique over the cohort.
10. **`driver.catalog.tsv` has no SV or cluster ID.** It is a gene list
    (PURPLE's reported drivers plus LINX disruptions), not an event
    table, so the report does not use it.
11. **LINX's artifact rules were undocumented here, and the README
    differs from the code.** They are now taken from the code (6.4).
12. **Two statements in `LINX_tables.md` were wrong** ("an SV with both
    ends outside genes has no row" in `breakend`; rows "for all
    transcripts LINX matched"). Corrected from the code (5.1, 5.2). The
    same wording was in the table descriptions of `scripts/postOA.R`,
    which go into `manifest.tsv`; corrected there too, so the
    `manifest.tsv` text changes at the next `postOA.sh` build (the
    current one still has the old text).
13. **LINX calls fusions of a gene with itself.** 96 LINX fusions in 12
    tumors join two breakends of the same gene (names like `GENE_GENE`),
    made by 70 DEL, 27 DUP and 4 INV events, none reported. They join
    two parts of one transcript, not two genes (section 3); the user,
    who treats them as cryptic splice products, decided to
    leave them out of the fusion sheet, named `fusions_real`, and of
    the gene fusion counts, and to keep them on the events sheets
    (section 3).
14. **The genes sheet was too long to read** (8,224 genes, most in one
    tumor). It now lists genes with 2 or more tumors in at least one
    count (2,435).

## 11. Results of this run

SV table of `postOA.sh` job 13376226 (Prism commit `180d01f`), 14
groups (9 with a matched normal, 5 tumor-only), report rebuilt on
2026-11-23:

- 12,216 PASS events; 4,180 high confidence (section 6.6 for the
  flags).
- 6,938 events with a `start_gene`, 6,068 with an `end_gene`, 5,782
  with a disrupted gene at either end, 1,172 in a LINX fusion (1,073 in
  a fusion of two different genes), 9 in a fusion LINX reports.
- Same-gene fusions (left out of `fusions_real` and the gene fusion
  counts): 96 fusions on 101 events in 12 tumors.
- `genes` 2,435 rows, `fusions_real` 1,066 rows.

## 12. Open points

- The 212 PASS BND and SGL events that pass LINX's copy-number and
  distance tests but are not LOW_VAF and not LINE (6.4) were not
  traced.
- The flag rule is the user's choice of 2026-11-23; changing it means
  editing `flag_*` in `svReport.R` and the `rule` sheet text together
  (the glossary points to the rule sheet).
- A nearest-gene column is not provided (decision 2026-11-23).
