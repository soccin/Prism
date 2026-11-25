# Prism post-pipeline methods

This file records what Prism does after the oncoanalyser pipeline has
finished. Each post-pipeline step has its own section, added when the
step is built. The definitions are taken from the code at the time of
writing, not from READMEs or summaries, so that a client or the methods
section of a paper can rely on them. When the code changes, change the
section with it.

What the pipeline itself does (oncoanalyser 3.0.0 and the Hartwig
tools) is in `ACTUAL_METHODS.md`. This file starts where that one ends:
the input of every step here is the finished `out/` folder.

## Project rule

Post-pipeline steps do not infer biology. They preserve what the tools
wrote:

- no record is removed;
- every value a tool wrote can be recovered from the output exactly
  (the SV table reshapes some fields for analysis, section 1.2; the
  check in section 1.4 rebuilds every source value from it);
- anything added is a documented extra column;
- a confidence tier, if one is ever added, is a column, never a filter;
- an error stops the run; nothing fails silently.

The SV report (section 2) has a confidence column and two user-decided
exceptions in its summary sheets only (same-gene LINX fusions left out
of the fusion sheet and gene fusion counts; genes in fewer than 2
tumors left out of the genes sheet); its events sheet keeps every PASS
event with every value it shows.

## Software

| Software | Version | Source of the version |
|---|---|---|
| nf-core/oncoanalyser | 3.0.0 (`v3.0.0-g7c74c87`), Nextflow 25.10.4 | `out/pipeline_info/nf_core_oncoanalyser_software_mqc_versions.yml` |
| ESVEE | 2.0.1 | same file |
| PURPLE | 4.4 | same file |
| LINX | 2.3.1 | same file; also `samples.tsv` column `linx_version` |
| R | 4.5.1 | `Rscript` on the `PATH` of the job |
| Bioconductor | 3.22: GenomicRanges 1.62.1, IRanges 2.44.0, S4Vectors 0.48.1 | same |
| CRAN | tidyverse 2.0.0 (dplyr 1.2.1, readr 2.2.0, tidyr 1.3.2, stringr 1.6.0, purrr 1.2.2), fs 2.1.0, nanoparquet 0.5.2, openxlsx 4.2.9 (`svReport.R`, `orangeSvReport.R`), jsonlite 2.0.0 (`orangeSvReport.R`) | same |
| Python | 3.10.14, pyarrow 25.0.1 | `python3` on the `PATH` of the job |
| pysam | 0.24.0 (its own htslib 1.23.1) | `python3` on the `PATH` of the `bamEvidence.sh` job |
| bcftools | 1.20 | `/admin/software/bcftools/bcftools-1.20/bin/bcftools` (`BCFTOOLS` in `popSvSites.sh` and `spotCheckSvCohort.R`); first line of the `popSvSites.sh` log |
| samtools | 1.20 (htslib 1.20) | `samtools` on the `PATH`; used by `spotCheckSvCohort.R` only |

The R and pyarrow versions were read on 2026-11-22, pysam, samtools
and openxlsx on 2026-11-23, from the `Rscript`, `python3` and `samtools` on the
`PATH` of the jobs. `postOA.sh` and `bamEvidence.sh` do not log them.
To record them for a run, use the same environment as the job and run
(a few seconds):

```
Rscript --vanilla -e "cat(R.version.string, '\n'); for (p in c('BiocVersion', 'GenomicRanges', 'IRanges', 'S4Vectors', 'tidyverse', 'dplyr', 'readr', 'tidyr', 'stringr', 'purrr', 'fs', 'nanoparquet')) cat(p, as.character(packageVersion(p)), '\n')"
python3 -c "import sys, pyarrow; print(sys.version.split()[0], 'pyarrow', pyarrow.__version__)"
python3 -c "import pysam, pysam.version; print('pysam', pysam.__version__, 'htslib', pysam.version.__htslib_version__)"
samtools --version
```

The Prism version of a run is the `PRISM=` line at the top of its log
(`git describe --always --dirty` of the Prism folder).

## 1. postOA.sh: cohort tables

```
./Prism/scripts/fetchTracks.sh [TRACK_DIR]           # once per site, login node
sbatch ./Prism/scripts/popSvSites.sh [TRACK_DIR]     # once per site, after fetchTracks.sh
./Prism/bamEvidence.sh [OUT_DIR]                     # per run: prints the sbatch command below
sbatch --array=0-N ./Prism/bamEvidence.sh [OUT_DIR]  # per run, N = number of BAMs - 1
sbatch ./Prism/postOA.sh [-t TRACK_DIR] [OUT_DIR]
```

Three inputs of step 2 are made before `postOA.sh` runs:

- `scripts/fetchTracks.sh` downloads the region tracks and the
  population SV call sets into `TRACK_DIR` (once per site; see
  Coordinates and region tracks, and Population SV call sets);
- `scripts/popSvSites.sh` reduces the population SV VCFs to one site
  table per source (once per site; see Population SV call sets);
- `bamEvidence.sh` counts the reads at every SV breakend of the cohort
  in every BAM of the run (once per run; see BAM evidence). It reads
  `res/prism/db/samples.tsv`, which step 1 writes. On a new run, run
  `postOA.sh` once (it stops at step 2: `--bam-dir not found`) or
  `Rscript Prism/scripts/postOA.R --out-dir OUT_DIR --res-dir
  res/prism/db` alone, then `bamEvidence.sh`, then `postOA.sh`.

Section 1.5 is a spot check of the computed columns, run by hand after
a build; `postOA.sh` does not run it.

`postOA.sh` runs four steps. All write to `res/prism/db/` in the folder
that holds `OUT_DIR`:

1. `scripts/postOA.R`: the stacked LINX tables (section 1.1);
2. `scripts/svCohort.R`: the SV event table as text, with its column
   dictionary (section 1.2);
3. `scripts/tsvToParquet.py`: the typed Parquet copy of that table
   (section 1.3);
4. `docs/scripts/checkSvCohort.R`: rebuilds every source value from the
   table and checks the Parquet copy (section 1.4).

The job stops at the first step that fails (`set -e`; each step exits
with a non-zero status on an error). The log `postOA.<jobid>.out` is
moved to `OUT_DIR/logs/other/` at the end of the job. A rerun replaces
every file.

Before the first step, `postOA.sh` stops if `OUT_DIR` does not exist,
if `TRACK_DIR/tracks.tsv` does not exist, if `OUT_DIR/pipeline_info/`
has no `params_*.json`, if the Hartwig resource folder named there
(`ref_data_hmf_data_path` of the last file in name order) does not
exist, if the Prism folder cannot be found, if `Rscript` is not on the
`PATH`, if the R package `nanoparquet` is not installed, or if
`python3` cannot import `pyarrow`.

### 1.1 Step 1, `scripts/postOA.R`: stacked LINX tables

**Groups.** One row of `samples.tsv` per folder of `OUT_DIR` that holds
`purple/*.purple.purity.tsv`, in order of the folder name:

| Column | Definition |
|---|---|
| `group_id` | the folder name |
| `sample_id` | the tumor: the purity file name without `.purple.purity.tsv` (the first file if there are several) |
| `normal_id` | the other sample folders in `bamtools/`, joined with `;`; empty if there is none |
| `run_mode` | `runMode` of the purity file: `TUMOR_GERMLINE` (matched normal) or `TUMOR` (tumor-only) |
| `linx_version` | the `version=` line of `linx/somatic_annotations/linx.version` |

The script warns (it does not stop) when `run_mode` and `normal_id`
disagree. Step 2 treats a group as matched when `normal_id` is not
empty.

**Tables.** For each LINX table, the per-sample files of all groups are
stacked into `data/linx_<mode>_<table>.tsv.gz`:

- somatic, every group: `svs`, `clusters`, `links`, `breakend`,
  `fusion`, `drivers`, `driver_catalog`, from
  `linx/somatic_annotations/<sample_id>.linx.<table>.tsv`;
- germline, matched groups: `svs`, `breakend`, `disruption`,
  `clusters`, `links`, `driver_catalog`, from
  `linx/germline_annotations/<sample_id>.linx.germline.<table>.tsv`.

The output header is `group_id`, `sample_id` and the header of the
files. Each data line is `group_id`, `sample_id` and the line as LINX
wrote it: the lines are not parsed, so no value is changed, reformatted
or quoted. LINX's `neoepitope` and `vis_*` files are not stacked. What
one row of each table is, and how the tables link, is in
`LINX_tables.md` and `ACTUAL_METHODS.md` section 4.

**Stops.** No group with a purity file; any LINX file missing (all
missing files are listed); a file with no header line; a header that
differs from the first file of the table; a line with a different
number of tab-separated fields than its header.

**Other files.** `samples.tsv`; `rows_per_sample.tsv` (the rows each
group adds to each table; empty where the table does not apply);
`manifest.tsv` (one row per file: tool, mode, number of groups and
rows, source path, description).

### 1.2 Step 2, `scripts/svCohort.R`: the SV event table

The SV event table is the core analysis table of the project. Output:
`res/prism/db/sv_cohort_filtered.tsv.gz` (text) and its column
dictionary `sv_cohort_filtered.columns.tsv`. Step 3 writes the same
table as `sv_cohort_filtered.parquet`, typed; that is the file to
analyze (section 1.3).

**The name.** "filtered" refers to the input, not to anything Prism
does. The input is PURPLE's somatic SV file `purple.sv.vcf.gz`, which
holds the calls ESVEE kept after its own filtering (FILTER `PASS` and
`PON`), plus PURPLE's inferred breakends (`INFERRED`); see
`ACTUAL_METHODS.md` section 3. Prism removes nothing: every record of
that file is in the table, whatever its FILTER value.

#### Inputs

| Input | Path |
|---|---|
| groups | `res/prism/db/samples.tsv` (step 1) |
| events | `out/<group_id>/purple/<sample_id>.purple.sv.vcf.gz` |
| ESVEE breakend table | `out/<group_id>/esvee/<sample_id>.esvee.breakend.tsv` |
| LINX tables | `res/prism/db/data/linx_somatic_{svs,clusters,links,breakend,fusion,drivers}.tsv.gz` (step 1) |
| cross-sample sets | `purple.sv.vcf.gz`, `purple.sv.germline.vcf.gz`, `esvee/<sample_id>.esvee.unfiltered.vcf.gz` of every group |
| region tracks | `TRACK_DIR/tracks.tsv` and the files it lists |
| Hartwig regions | `<hmf>/dna/sv/sv_prep_blacklist.38.bed`, `<hmf>/common/unmap_regions.38.tsv` |
| driver gene panel | `<hmf>/common/DriverGenePanel.38.tsv` (the panel LINX reads; columns `gene`, `likelihoodType`) |
| BAM read counts | `res/prism/db/bam_evidence/<sample>.bam_evidence.tsv.gz`, one file per tumor and per normal of `samples.tsv` (`bamEvidence.sh`; see BAM evidence) |
| population SV sites | `TRACK_DIR/<source>/sites.tsv.gz` for `gnomad_sv`, `kgp_sv`, `dbvar_nstd166`, `dbvar_nstd200` (`scripts/popSvSites.sh`; see Population SV call sets) |

Every file is read as text (`col_character`, no quote handling, no
value turned into NA), so values reach the table exactly as written.
The one exception is the VCF missing value `.`: in QUAL, FILTER, INFO
and FORMAT it is stored as an empty cell. The BAM read counts and the
population sites are Prism's own intermediate files: their numbers are
read as numbers, `NA` (BAM counts) and `.` (sites) as missing.

#### Events (rows)

1. Each record of `purple.sv.vcf.gz` is one breakend. Its orientation
   is derived from ALT: 1 if ALT starts with a base (`A`, `C`, `G`,
   `T`, `N`, either case), else -1. This is the convention of ESVEE and
   LINX; the ESVEE check below confirms it for every record.
2. Two records form one event when each one's INFO `MATEID` is the
   other's ID. A `MATEID` that is not reciprocal stops the script. A
   record whose `MATEID` names an ID that is not in the file is treated
   as a one-breakend event.
3. The start of a two-breakend event is the record that comes first in
   the VCF; the other is the end.
4. A record without a mate (a single breakend, or a breakend PURPLE
   inferred from copy number, FILTER `INFERRED`) is a one-breakend
   event: it is the start, and every `end_*` column is empty, except
   `end_linx_bnd_*` where LINX wrote end rows for it (see Joins).
5. The script stops unless 2 x (two-breakend events) + (one-breakend
   events) = the number of records, so every record is in exactly one
   row.

Events are listed group by group (order of `samples.tsv`), and within a
group in VCF order of their start record.

#### Layout

Column order:

1. the lead block (user decision, table below): the event, the genes
   LINX attaches to each breakend (genes containing it or starting up
   to 100 kb after it; `SV_REPORT.md` section 5) with their role and
   LINX's annotation, fusions, cluster drivers and category, copy
   number and support;
2. the other key columns (table further below);
3. event-level columns: fields with the same value on both breakends,
   stored once, without a `start_`/`end_` prefix (see Unpacking, rule
   5);
4. every `start_*` column, then every `end_*` column (the same set of
   columns for each breakend);
5. the per-event LINX columns (`linx_*`);
6. the `xs_pair_*` columns, then the per-event population SV columns.

A column appears once: the lead block and key columns are left out of
the groups after them.

Lead block, in this order (53 columns):

| # | Columns |
|---|---|
| 1-5 | `group_id`, `sample_id`, `event_key`, `run_mode`, `filter` |
| 6-7 | `linx_type` (empty for PON events, which LINX does not analyze), `info_SVTYPE` (every event; BND for INFERRED, which LINX types INF) |
| 8-14 | `start_chrom`, `start_pos`, `start_orient`, `end_chrom`, `end_pos`, `end_orient`, `length` |
| 15-18 | `start_info_CIPOS_lower`, `_upper`, `end_info_CIPOS_lower`, `_upper` |
| 19-28 | `start_linx_bnd_gene`, `start_gene_role`, `start_linx_bnd_` `canonical`, `regionType`, `codingType`, `disruptive`, `exonUp`, `exonDown`, `undisruptedCopyNumber`, `reportedStatus` |
| 29-38 | the same 10 for `end_` |
| 39-43 | `fusion_gene_5p`, `fusion_gene_3p`, `linx_fusion_phased`, `linx_fusion_reported`, `linx_fusion_reportedType` |
| 44-48 | `linx_cluster_driver_gene`, `linx_cluster_driver_eventType`, `linx_cluster_category`, `linx_cluster_resolvedType`, `linx_cluster_clusterCount` |
| 49-53 | `start_info_PURPLE_CN_CHANGE`, `end_info_PURPLE_CN_CHANGE`, `info_PURPLE_JCN`, `start_info_PURPLE_AF`, `tumor_VF` |

The other `start_linx_bnd_*`, `end_linx_bnd_*` and `linx_fusion_*`
list columns stay in their groups; their elements line up by position
with the lists in the lead block.

Within a breakend, `info_*`, `tumor_*`, `esvee_bnd_*` and `normal_*`
columns appear in the order their keys are first seen, group by group,
and a split pair takes the place of its field. So the order can differ
between runs. Select columns by name.

Key columns (`group_id`, `sample_id`, `run_mode`, `event_key` and
`length` are in the lead block; the rest follow it in this order:
`normal_id`, `n_breakends`, `in_linx`, `linx_start_is_vcf_start`,
`same_chrom`, `start_coords`, `end_coords`, `start_id`, `end_id`):

| Column | Definition |
|---|---|
| `group_id`, `sample_id`, `normal_id`, `run_mode` | from `samples.tsv` |
| `event_key` | `sample_id:start_id` |
| `n_breakends` | 2 for a two-breakend event, 1 otherwise |
| `in_linx` | TRUE if LINX `svs.tsv` has a row for the event (see Joins) |
| `linx_start_is_vcf_start` | `linx_coordsStart` = `start_coords`; empty if not in LINX |
| `same_chrom` | both breakends on one chromosome; empty for one-breakend events |
| `length` | \|end POS - start POS\| when `same_chrom` is TRUE; else empty |
| `start_coords`, `end_coords` | `chromosome:position:orientation` |
| `start_id`, `end_id` | VCF IDs of the start and end records |

Per breakend (`start_` or `end_` in front):

| Columns | Source |
|---|---|
| `chrom`, `pos`, `ref`, `alt` | the VCF fields |
| `orient` | derived from ALT (above) |
| `vcf_order` | record number in the VCF |
| `info_<KEY>` | one column per INFO key seen in any record of any group (unless the key is event-level or a dropped repeat); the value after `=`, or `true` for a flag; empty where the record has no such key |
| `tumor_<KEY>`, `normal_<KEY>` | one column per FORMAT key, for the tumor sample (column `sample_id`) and the normal (the other sample column, if any); empty where the record has no value for the key |
| `esvee_bnd_<column>` | the columns of ESVEE's `breakend.tsv` that no VCF field repeats (see Joins and Unpacking) |
| `linx_bnd_<column>`, `linx_bnd_n` | LINX `breakend.tsv` rows at this breakend (every column but `svId`, including LINX's `coords`), one list element per row (see Joins and Lists) |
| `gene_role` | one element per element of `linx_bnd_gene`, same order: the gene's `likelihoodType` (ONCO or TSG) in the driver gene panel; empty element for a gene not on the panel (see Joins) |
| `region_*` | region tracks (see Region tracks) |
| `xs_*`, `own_germline_*` | the same breakend in other files (see Cross-sample columns) |
| `bam_tumor_*`, `bam_normal_*` | reads at the breakend recounted in the group's own tumor and normal BAM (see BAM evidence) |
| `bam_other_normals_*`, `bam_other_tumors_*` | the same counts summarized over the normal and tumor BAMs of the other groups (see BAM evidence) |
| `<source>_w100_n`, `<source>_w100_max_af` | sites of a population SV call set with a breakpoint within 100 bases (see Population SV call sets) |

Per event: `linx_<column>` (LINX `svs.tsv`), `linx_cluster_<column>`
(`clusters.tsv`), `linx_fusion_*`, `fusion_gene_5p`,
`fusion_gene_3p`, `linx_link_*`, `linx_cluster_driver_*` (lists),
`xs_pair_*`, and the population SV call set matches of the event,
`<source>_bp100_*` (by both breakpoints) and `<source>_ro50_*` (by
reciprocal overlap) (see Population SV call sets). `<source>` is
`gnomad_sv`, `kgp_sv`, `dbvar_nstd166` or `dbvar_nstd200`.
`gene_role` and `fusion_gene_*` do not carry the
`linx_bnd_` or `linx_fusion_` prefix, which is kept for columns of the
LINX files (`checkSvCohort.R` rebuilds the LINX rows from them).

In the TSV, logical columns Prism adds are written `TRUE` or `FALSE`;
values copied from a tool keep the tool's spelling (LINX writes
`true`/`false`). In the Parquet copy both are boolean.

#### Joins

- **VCF to ESVEE.** VCF ID = ESVEE `Id`. Only ESVEE rows whose `Id` is
  in the PURPLE file are used. For every record except FILTER
  `INFERRED`, the ESVEE row must exist and have the same chromosome,
  position and orientation as the VCF record; otherwise the script
  stops. Inferred breakends are PURPLE's and have no ESVEE row; their
  `esvee_bnd_*` cells are empty. Most ESVEE columns repeat a VCF field
  and are dropped after a check (Unpacking, rule 1).
- **Event to LINX `svs.tsv`.** `vcfIdStart` = `start_id`, within the
  same `group_id` and `sample_id`. For each match, `vcfIdEnd` must be
  `end_id` (empty for a one-breakend event) or equal `start_id`; LINX
  repeats the start ID in `vcfIdEnd` for some SVs whose two breakends
  share a position (`ACTUAL_METHODS.md` section 4). Every LINX `svs`
  row must match an event. LINX analyzes only FILTER `PASS` and
  `INFERRED`, so PON events have `in_linx` FALSE and empty LINX
  columns.
- **LINX start and end.** For some SVs whose two breakends share a
  position, LINX's start breakend is the VCF end
  (`linx_start_is_vcf_start` FALSE; LINX's `coordsStart`, `vcfIdStart`
  and other start columns are kept as LINX wrote them). If LINX's start
  differs from the VCF start for an SV whose breakends are at different
  positions, the script stops.
- **Clusters.** `linx_clusterId` = `clusters.clusterId`.
- **LINX `breakend.tsv`.** A breakend row (one per transcript) goes to
  the event with the same `svId`, and to one of its two ends by LINX's
  `isStart`. A row with `isStart` `true` goes to the start
  (`start_linx_bnd_*`) and a row with `false` to the end
  (`end_linx_bnd_*`), except when `linx_start_is_vcf_start` is FALSE,
  where the two are swapped. LINX's `coords` of each row is kept
  (`start_linx_bnd_coords`, `end_linx_bnd_coords`). Rows are not placed
  by coordinates, for two reasons. An SV can have both breakends at the
  same chromosome, position and orientation, so a match by coordinates
  would place its rows on both ends. And for some single breakends
  LINX writes `isStart` `false` rows at a second location that is not a
  VCF record (empty `vcfId`). Those rows fill `end_linx_bnd_*` of a
  one-breakend event whose `end_*` VCF columns are empty. The script
  stops if a breakend row's `svId` matches no LINX event, if `isStart`
  is not `true` or `false`, or if the placed rows do not number exactly
  the rows of the table. So each row is placed exactly once.
- **Fusions.** Each fusion's `fivePrimeBreakendId` and
  `threePrimeBreakendId` are looked up in `breakend.id` to get the
  SV(s). A fusion is listed on the SV of its 5' breakend with
  `linx_fusion_side` = `5p` and on the SV of its 3' breakend with `3p`.
  When both breakends are on one SV, that SV lists the fusion twice
  (`5p`, `3p`). A fusion breakend ID without a `breakend.tsv` row
  stops the script.
- **Fusion genes.** `fusion_gene_5p` and `fusion_gene_3p` are the
  `breakend.tsv` `gene` of `fivePrimeBreakendId` and
  `threePrimeBreakendId`, one element per fusion in the order of the
  `linx_fusion_*` lists. They are not split from `linx_fusion_name`,
  because the name does not always hold the breakend genes. The script
  checks every fusion: the name must be `<5' gene>_<3' gene>`, or, for
  an enhancer partner on the 5' side, `@<region>_<3' gene>`; anything
  else stops it. LINX names enhancer partners by region (`@IGH`,
  `@IGK`, `@BCL11B`) while the breakend's gene is another: in this
  cohort 12 of 1,171 fusions, all `ENHANCER_PROMISCUOUS`, e.g.
  `@IGH_PSD3` with 5' breakend gene IGHA1, `@BCL11B_LRP12` with RPL3P4.
- **Gene roles.** Each `breakend.tsv` row's `gene` is looked up in the
  driver gene panel (`DriverGenePanel.38.tsv`; 514 genes: 207 ONCO, 307
  TSG). The roles are collapsed from the same rows in the same order as
  `{side}_linx_bnd_*`, so `{side}_gene_role` has as many elements as
  `{side}_linx_bnd_n` (checked on every row, then checked again element
  by element against the panel). The script stops if the panel is
  missing, lists a gene twice, or has a `likelihoodType` other than
  ONCO or TSG. The roles are Hartwig's pan-cancer calls, mostly from
  solid tumors: JAK1 and ERBB4 are TSG there, while JAK1 is activated
  in some T-cell lymphomas. The role says nothing about the effect of
  this event.
- **Links.** A link is listed on its `lowerSvId` and, if different, on
  its `upperSvId`.
- **Cluster drivers.** `drivers.tsv` rows are listed on every event of
  their cluster (`clusterId` = `linx_clusterId`). Rows with
  `clusterId` -1 (no cluster) are on no event.

#### Unpacking

The tools write some fields as packed values: two numbers in one VCF
field, the same value on both records of an SV, the same value in the
VCF and in ESVEE's table. For analysis, a cell holds one value, except
where the source holds a list of any length. Five rules, applied in
this order, reshape the table. Each rule checks its assumption on every
row first and stops the run, naming the column and an example value, if
a row breaks it. Section 1.4 undoes every rule and compares the result
with the source files.

**Rule 1, repeated columns are dropped.** A column is dropped only after
every value it holds is shown to equal the kept column (an empty value
in the dropped column is allowed, for example ESVEE for an `INFERRED`
record). Per breakend; "other" means the other breakend of the event.

| Dropped | Equal to | Note |
|---|---|---|
| `esvee_bnd_AssemblyId` | `info_ASMID` | |
| `esvee_bnd_AssemblyInfo` | `info_ASMINFO` | |
| `esvee_bnd_Type` | `info_SVTYPE` | |
| `esvee_bnd_Chromosome`, `_Position`, `_Orientation` | `chrom`, `pos`, `orient` | |
| `info_BEOR` | `orient` | |
| `esvee_bnd_MateId`, `info_MATEID` | the other breakend's `id` | |
| `esvee_bnd_MateChr`, `_MatePos`, `_MateOrient` | the other breakend's `chrom`, `pos`, `orient` | |
| `esvee_bnd_ConfidenceInterval` | `info_CIPOS` | |
| `esvee_bnd_InexactOffset` | `info_IHOMPOS` | |
| `esvee_bnd_Homology` | `info_HOMSEQ` | |
| `esvee_bnd_Qual` | `qual` | |
| `esvee_bnd_SplitFragments`, `_DiscFragments` | `info_SF`, `info_DF` | |
| `esvee_bnd_SequenceLength` | `info_ASMLEN` | |
| `esvee_bnd_SegmentIndex`, `_SequenceIndex` | `info_ASMSEG`, `info_BEAPOS` | |
| `esvee_bnd_AlignedBases`, `_MapQual`, `_Score`, `_AdjAlignedBases` | `info_SEGALEN`, `info_SEGMAPQ`, `info_SEGSCO`, `info_SEGRL` | |
| `esvee_bnd_FacingBreakendIds` | `info_ASMLNKS` | ESVEE separates with `;`, the VCF with `,` |
| `esvee_bnd_AltAlignments` | `info_INSALN` | |
| `esvee_bnd_UniqueFragPos` | `info_UFP` | |
| `esvee_bnd_ThreePrimeRange` | `info_TPR` | ESVEE separates with `;`, the VCF with `,` |
| `esvee_bnd_StrandBias` | `tumor_SB`, then `normal_SB`, joined with `;` | ESVEE writes SB rounded to two decimals (`%.2f` in ESVEE's `BreakendWriter`); equal within 0.005 |
| `vcf_normal_sample` | `normal_id` | the VCF's second sample column is the normal of `samples.tsv` |

ESVEE `AvgFragmentLength` is not a repeat of INFO `AVGLEN` (the values
differ) and is kept.

**Rule 2, PURPLE's per-breakend pairs become one value per breakend.**
PURPLE writes `PURPLE_AF`, `PURPLE_CN` and `PURPLE_CN_CHANGE` as two
values: [start, end] on the start record and [end, start] on the end
record (`SomaticSvCache` at `purple-v4.4`). The rule checks that, for
every two-breakend event, both records hold exactly two values and the
second value of each is the first value of the other (or both are
empty); for a one-breakend event, that the start holds at most one
value. It then keeps each record's first value, its own:
`start_info_PURPLE_AF` is the start's value, `end_info_PURPLE_AF` the
end's.

**Rule 3, fixed pairs are split into two columns.** Fields that always
hold exactly two values (checked on every row) become two columns,
named `<field>_<part>`:

| Field | Parts | Meaning |
|---|---|---|
| `info_CIPOS` | `lower`, `upper` | confidence interval around POS: lower and upper offset |
| `info_IHOMPOS` | `lower`, `upper` | inexact homology: lower and upper offset |
| `info_TPR` | `five_prime_span`, `three_prime_span` | spread (max - min) of the 5' and of the 3' positions of the unpaired reads at the breakend, -1 if there are none (ESVEE `Breakend.readOrientationRange`); spans in bases, not positions |
| `tumor_AD`, `normal_AD` | `ref`, `alt` | reads supporting the reference and the alt allele |

**Rule 4, lists that always hold one value become scalars.** These
fields are lists by definition but hold one value in every row; the
rule checks that no value contains the separator and keeps the column
as one value. A row with more than one value stops the run.

| Field | Separator |
|---|---|
| `filter` | `;` |
| `info_ASMID`, `info_ASMSEG`, `info_BEAOR`, `info_BEAPOS`, `info_SEGID` | `,` |
| `region_encode_blacklist_reason` | `\|` |

**Rule 5, fields equal on both breakends are event-level.** For each
field below, the rule checks that every two-breakend event has the same
value on both records (or both empty) and that no one-breakend event has
an end value. The field then becomes one column without a prefix,
holding the start's value. A field that no record has gets no column.

`qual`, `filter`; INFO `SVTYPE`, `PURPLE_JCN`, `SVID`, `ASMID`,
`ASMINFO`, `ASMLEN`, `AVGLEN`, `SF`, `DF`, `VF`, `MLR`, `INSALN`,
`INSRMP`, `INSRMRC`, `INSRMRT`, `PON_COUNT`, `TINC_RECOVERED`; FORMAT
`GT`, `SF`, `DF`, `VF` of the tumor and the normal; ESVEE
`PhaseGroupId`, `PhaseSetId`, `Length`, `RefSplitFragments`,
`RefDiscFragments`, `SegmentCount`, `IncompleteFragments`,
`AvgFragmentLength`, `NonPrimaryFragments`. (The 2-group test has 35 of
these 36; it has no `info_TINC_RECOVERED`.)

#### Lists

A list column holds 1 to N values in each non-empty cell. In the
Parquet copy it is an array; in the TSV the elements are joined with
the column's separator (`list_sep` in the dictionary):

| Columns | Separator | Elements |
|---|---|---|
| `start_linx_bnd_*`, `end_linx_bnd_*` | `\|` | one per LINX `breakend.tsv` row placed on the breakend |
| `linx_fusion_*`, `linx_link_*`, `linx_cluster_driver_*` | `\|` | one per fusion, link or driver row of the event |
| `region_rmsk_name`, `_class`, `_family` | `\|` | distinct labels of the repeats at the position |
| `*_group_list_w50` | `\|` | distinct group IDs, sorted |
| `<source>_bp100_ids`, `<source>_ro50_ids` | `\|` | IDs of the matching population sites, each once; a set, the order carries no meaning |
| `info_*` with VCF Number other than 0 or 1 (for example `ASMLNKS`) | `,` | as written in the VCF |
| `linx_clusterReason`, `linx_lineTypeStart`, `linx_lineTypeEnd`, `linx_geneStart`, `linx_geneEnd` | `;` | as LINX wrote them in `svs.tsv` |

Lists built from source rows (the first two rows of the table):

- The columns of one source table are aligned by position: element k
  of every `start_linx_bnd_*` cell of a row comes from the same LINX
  breakend row, and the same holds for `end_linx_bnd_*`,
  `linx_fusion_*`, `linx_link_*` and `linx_cluster_driver_*`.
- Order: `breakend.tsv`, `drivers.tsv`: file order. Fusions: first the
  fusions where the SV is the 5' side, then those where it is the 3'
  side, each in file order. Links: first the links where the SV is
  `lowerSvId`, then those where it is `upperSvId`, each in file order.
- `<prefix>n` (`start_linx_bnd_n`, `end_linx_bnd_n`, `linx_fusion_n`,
  `linx_link_n`, `linx_cluster_driver_n`, the dictionary's
  `list_count_column`) holds the number of source rows, so every list
  of the family has exactly n elements. A fusion listed twice on one SV
  counts twice.
- An empty element is an empty source cell (null in Parquet). In the
  TSV an empty cell cannot be told from a list of one empty element;
  the count column decides.
- The script stops if any value to be joined contains `|`.

The RepeatMasker labels are sets, not aligned lists: the distinct
labels are taken in file order, separately in each column, so element k
of `region_rmsk_name` need not belong to element k of
`region_rmsk_class`.

The population site ID lists are not checked for `|` (the check above
is in the LINX list code only). None of the four site tables contains
a `|` (checked on 2026-11-23).

**Empty cells.** In the TSV an empty cell is a missing value; in the
Parquet copy it is null. A cell is empty when: the column is an
`end_*` column of a one-breakend event (except LINX end rows, see
Joins); the record has no such INFO or FORMAT key, or the VCF has `.`;
the group is tumor-only (`normal_*`, `own_germline_*`, `bam_normal_*`);
the record is `INFERRED` (`esvee_bnd_*`); the event is not in LINX
(`linx_*`); there is no source row for a list; there is no interval or
group to list; a BAM measure is not defined for the breakend (no mate,
no reads; see BAM evidence); or a population rule does not apply to the
event (`<source>_bp100_*` for one-breakend events, `<source>_ro50_*`
for events other than DEL, DUP and INV) or no matching site has an AF
(`*_max_af`). A count is 0, not empty, where the rule applies and
nothing matches.

#### Types

The dictionary gives each column one type: `string`, `integer`,
`double`, `boolean`, or `list<string|integer|double|boolean>`. The
type is set as follows, in this order:

1. Prism's own columns have a fixed type (for example `pos`, `orient`,
   `length`, `n_breakends`: integer; `qual`: double; `in_linx`,
   `same_chrom`: boolean; IDs and coordinates: string).
2. LINX `svs.tsv` fields that hold `;`-joined items (above):
   `list<string>`.
3. List families (`linx_bnd_*`, `linx_fusion_*`, `linx_link_*`,
   `linx_cluster_driver_*`): the count column is integer; the element
   type of each list is taken from all of its values, except VCF ID
   columns (`vcfId`, `*VcfId`), which are string.
4. VCF fields follow the header's `Number` and `Type`: `Number=0` is
   boolean; `Number=1`, a split part, `PURPLE_AF`/`CN`/`CN_CHANGE` and
   the rule 4 fields are one value; any other `Number` is a list on `,`.
   `Type` maps Integer to integer, Float to double, String and
   Character to string, Flag to boolean. One override: ESVEE declares
   `SEGID` Integer but writes `<assembly>_<segment>` (for example
   `32022_0`), so it is string. The script stops if a key's `Number` or
   `Type` differs between groups, or if a VCF column has no header
   definition.
5. Region columns: in-track flags boolean; `region_centromere_distance`
   integer; RepeatMasker labels `list<string>`;
   `region_encode_blacklist_reason` string. `xs_*`, `xs_pair_*` counts
   and `own_germline_*`: integer; group lists `list<string>`. BAM
   columns: counts integer (`depth`, `clip_n`, `split_n`, `disc_n`, and
   `n`, `*_samples`, `*_reads` of the other-sample summaries), means and
   fractions double. Population columns: `_n` and `_n_pass` integer,
   `_max_af` double, `_ids` `list<string>`.
6. Any other column (ESVEE and LINX scalars): taken from all of its
   values: boolean if all are `true`/`false`/`TRUE`/`FALSE`, integer if
   all are whole numbers, double if all are numbers, else string. A
   column with no values is string. Because this rule reads the values,
   such a column can get a different type in another cohort.

Step 3 converts every cell to its type and stops on any value that does
not fit, so a type in the dictionary is a checked statement about every
cell.

#### Coordinates and region tracks

VCF positions are 1-based. Every interval is turned into a 1-based
interval with both ends included before it is compared with a breakend
position:

| File | Format | Interval used |
|---|---|---|
| GIAB, ENCODE (`.bed.gz`) | BED: 0-based start, end excluded | [start + 1, end] |
| UCSC `rmsk.txt.gz` | `genoStart` 0-based, `genoEnd` excluded (columns 6 to 8) | [genoStart + 1, genoEnd] |
| UCSC `centromeres.txt.gz` | `chromStart` 0-based, `chromEnd` excluded (columns 2 to 4) | [chromStart + 1, chromEnd] |
| Hartwig `sv_prep_blacklist.38.bed` | header `Chromosome`, `PosStart`, `PosEnd`; `.bed`, so 0-based | [PosStart + 1, PosEnd] |
| Hartwig `unmap_regions.38.tsv` | header `Chromosome`, `PosStart`, `PosEnd`; not `.bed`, so 1-based | [PosStart, PosEnd] |

The Hartwig files are read the way Hartwig's own region loader reads
them (`ChrBaseRegion.loadChrBaseRegions` in hmf-common, tag
`esvee-v2.0.1`): it adds 1 to the start of any file whose name ends in
`.bed` or `.bed.gz`, header or not, and uses other files as written.

A breakend is in a track when its chromosome is the interval's and its
position lies in the interval. The point is a single base; no window
is added (some GIAB files already include flanks, named `slop` in the
file name).

| Column | Value |
|---|---|
| `region_<track>` (17 GIAB tracks, `encode_blacklist`) | TRUE or FALSE |
| `region_encode_blacklist_reason` | the ENCODE label (column 4) at the position; one label (rule 4) |
| `region_rmsk_name`, `_class`, `_family` | RepeatMasker `repName`, `repClass`, `repFamily` (columns 11 to 13) of the repeats at the position (lists) |
| `region_hmf_sv_prep_blacklist`, `region_hmf_unmap` | TRUE or FALSE (see below) |
| `region_centromere_distance` | see below |

`region_hmf_sv_prep_blacklist` records whether the breakend is in
Hartwig's `sv_prep_blacklist.38.bed`, not what any tool did with it. In
this run no task reads that file: it is in no `command.sh` under
`out/logs/`, the oncoanalyser 3.0.0 ESVEE module passes only
`-unmap_regions`, and ESVEE 2.0.1's `PrepConfig` has no blacklist
option. `region_hmf_unmap` is the file ESVEE was given.

`region_centromere_distance` is computed with
`IRanges::distanceToNearest` against the UCSC centromere model
intervals (109 rows: several per chromosome, 24 chromosomes). It is the
distance to the nearest model interval on the same chromosome, on
either side; arms are not distinguished. For an interval [s, e] and a
position p: 0 if s <= p <= e; p - e - 1 if p > e; s - p - 1 if p < s.
So a position directly next to an interval is also 0, and a position
in a gap between two model intervals of one centromere gets the
distance to the nearer one. The cell is empty when the chromosome has
no model interval (for example chrM or an unplaced contig).

The 20 region tracks, as downloaded by `scripts/fetchTracks.sh` on
2026-11-21 (`TRACK_DIR/tracks.tsv` holds the same list with full URLs,
and also the four population SV call sets, listed under Population SV
call sets). `md5 check` is the comparison with GIAB's published md5
list (`GRCh38-genome-stratifications-md5s.txt`); the other sources of
region tracks publish none, so their md5 identifies the file used.

GIAB = `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all`,
UCSC = `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database`.

| Track | Source, version | File | md5 | md5 check |
|---|---|---|---|---|
| `giab_alldifficult` | GIAB v3.6, union of all difficult regions | GIAB/Union/GRCh38_alldifficultregions.bed.gz | 55a6b84934fcea76eebf38fd94503f9b | ok |
| `giab_lowmap_segdup` | GIAB v3.6, union of low mappability and segmental duplications | GIAB/Union/GRCh38_alllowmapandsegdupregions.bed.gz | b1b8a7e364aa8bc6b470b3b549e126d6 | ok |
| `giab_lowmappability` | GIAB v3.6, low mappability | GIAB/Mappability/GRCh38_lowmappabilityall.bed.gz | 8a78766bb99abf1d25a5b307d49a0324 | ok |
| `giab_segdups` | GIAB v3.6, segmental duplications | GIAB/SegmentalDuplications/GRCh38_segdups.bed.gz | 72ec3f7db18ef9318c64b1df271703fe | ok |
| `giab_tandem_homopolymer` | GIAB v3.6, tandem repeats and homopolymers, 5 bp flanks | GIAB/LowComplexity/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz | 05a302a437da3f3514106d3eedd62696 | ok |
| `giab_satellites` | GIAB v3.6, satellites, 5 bp flanks | GIAB/LowComplexity/GRCh38_satellites_slop5.bed.gz | 63da9414e5edad77a2b29ae93a47ab20 | ok |
| `giab_gc_extreme` | GIAB v3.6, GC below 25% or above 65%, 50 bp flanks | GIAB/GCcontent/GRCh38_gclt25orgt65_slop50.bed.gz | cd5fad19a02695f3e402c94a1874cb1c | ok |
| `giab_gaps` | GIAB v3.6, reference gaps, 15 kb flanks | GIAB/OtherDifficult/GRCh38_gaps_slop15kb.bed.gz | 194125f93992e1167d1b534fd23171c6 | ok |
| `giab_other_difficult` | GIAB v3.6, union of other difficult regions | GIAB/OtherDifficult/GRCh38_allOtherDifficultregions.bed.gz | 2303d01c4e32914dc424c90ec1e595a5 | ok |
| `giab_false_dup_incorrect` | GIAB v3.6, false duplications, incorrect copy | GIAB/OtherDifficult/GRCh38_false_duplications_incorrect_copy.bed.gz | 24112e2c656aa4b610a8e7d0f9012979 | ok |
| `giab_false_dup_correct` | GIAB v3.6, false duplications, correct copy | GIAB/OtherDifficult/GRCh38_false_duplications_correct_copy.bed.gz | 89e87651a0239b892170d7aad942f5fd | ok |
| `giab_collapsed_dup` | GIAB v3.6, collapsed duplication false-positive regions | GIAB/OtherDifficult/GRCh38_collapsed_duplication_FP_regions.bed.gz | 483cd58b013ff0ffb136e0700e5ee1bb | ok |
| `giab_population_cnv_fp` | GIAB v3.6, population CNV false-positive regions | GIAB/OtherDifficult/GRCh38_population_CNV_FP_regions.bed.gz | d74487e99d9f8a7cd79e14fe2a1987d2 | ok |
| `giab_l1h` | GIAB v3.6, L1H elements over 500 bp | GIAB/OtherDifficult/GRCh38_L1H_gt500.bed.gz | a0de22565c7a53834e1e267f2620d59a | ok |
| `giab_vdj` | GIAB v3.6, VDJ regions | GIAB/OtherDifficult/GRCh38_VDJ.bed.gz | bfe48ba5586f624972b8c243d21d97a0 | ok |
| `giab_mhc` | GIAB v3.6, MHC | GIAB/OtherDifficult/GRCh38_MHC.bed.gz | 47675293278e690d7ffce11fd6dafc2d | ok |
| `giab_kir` | GIAB v3.6, KIR | GIAB/OtherDifficult/GRCh38_KIR.bed.gz | a88877e50dad3a7bb567549e611f04cf | ok |
| `encode_blacklist` | ENCODE blacklist v2 (Boyle-Lab, `master` branch) | https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz | 83fe6bf8187a64dee8079b80f75ba289 | none |
| `ucsc_rmsk` | UCSC hg38 RepeatMasker table (no version; date and md5) | UCSC/rmsk.txt.gz | b2e108b535550ba9e3cf83c77417380f | none |
| `ucsc_centromeres` | UCSC hg38 centromere model table (no version; date and md5) | UCSC/centromeres.txt.gz | 94cb08b0b00932436cf4d5704257e224 | none |

The two Hartwig files are not downloaded; they come from the run's own
resource folder (`hmf_pipeline_resources.38_v3.0.0--8` in this
project), found through `ref_data_hmf_data_path` in the run's
`params_*.json`.

#### Cross-sample columns

These columns say whether the same breakend, or the same pair of
breakends, is in the calls of other groups. They are counts; they do
not mark or remove anything.

**Sets.** Three sets of breakends, read from every group:

| Set | Records |
|---|---|
| `somatic` | `purple.sv.vcf.gz`, every record except FILTER `INFERRED` (so PASS and PON) |
| `germline` | `purple.sv.germline.vcf.gz`, every record, of matched groups only (`normal_id` not empty) |
| `unfiltered` | `esvee.unfiltered.vcf.gz`, every record, any FILTER |

**Points.** Every breakend of every event: the start, and the end of
two-breakend events. This includes PON and `INFERRED` breakends.

**Hit.** A set record is a hit for a point when it is on the same
chromosome, has the same orientation (derived from ALT in the same
way), and |set POS - point POS| <= 50. "Exact" means the difference is
0. Positions only are compared: no confidence interval, homology or
insert sequence is used.

Per breakend (`start_` or `end_` in front), for each set S:

| Column | Definition |
|---|---|
| `xs_<S>_groups_exact` | number of other groups with an exact hit in their S records |
| `xs_<S>_groups_w50` | number of other groups with a hit within 50 bases |
| `xs_<S>_group_list_w50` | those groups (`group_id`), sorted (a list); empty if none |
| `own_germline_exact`, `own_germline_w50` | number of records (not groups) in the group's own `purple.sv.germline.vcf.gz` with an exact hit, or a hit within 50 bases |

"Other groups" excludes the event's own group in every `xs_*` column.
The counts are 0 where there is no hit. `own_germline_*` is empty for
tumor-only groups (no `normal_id` in `samples.tsv`) and 0 or more for
every matched group, including one whose germline file has no records.

Per two-breakend event, for each set S:

| Column | Definition |
|---|---|
| `xs_pair_<S>_groups_w50` | number of other groups that have two S records r1 and r2 with r1 a hit for the event's start, r2 a hit for the event's end (each within 50 bases, same orientation), and r1's `MATEID` = r2's ID |
| `xs_pair_<S>_group_list_w50` | those groups, sorted (a list); empty if none |

Both ends must match, and the two matching records must be one SV in
the other group. `xs_pair_*` is empty for one-breakend events and 0 for
two-breakend events without such a group.

#### BAM evidence

These columns count the reads at each breakend directly in the BAM
files, separately from the SV caller, in the group's own tumor and
normal and in the samples of the other groups. They are counts; they
do not mark or remove anything.

**Recount, `bamEvidence.sh` and `scripts/bamEvidence.py`.** One SLURM
array task per BAM: every `OUT_DIR/<group>/alignments/<sample>/<sample>.redux.bam`
(the REDUX output, duplicates marked), in sorted path order; in this
project 23 BAMs (14 tumors, 9 normals). Each task counts its BAM at
every breakend of the cohort: every record of every group's
`purple.sv.vcf.gz`, all FILTER values (32,965 records), not only the
records of its own group. Output:
`res/prism/db/bam_evidence/<sample>.bam_evidence.tsv.gz`, one row per
record: `group_id`, `id` (VCF ID) and the nine measures below, `NA`
where a measure is not defined, fractions and means to 6 significant
digits.

A breakend has its POS, its CIPOS (0,0 when absent), its orientation (1
if ALT starts with REF, -1 if ALT ends with REF) and its mate, the
record its `MATEID` names (none for single and inferred breakends).
Orientation 1 keeps the reference up to POS: supporting reads align up
to POS and continue elsewhere. Orientation -1 keeps the reference from
POS on.

Reads: unmapped, secondary, duplicate and QC-fail records are skipped.
Primary records give the depth, the clip fraction and the read pairs;
primary and supplementary records give the clipped reads. `clip_n`,
`split_n` and `disc_n` count distinct read names (fragments); `depth`
and the fractions count primary records. No mapping-quality threshold
is applied; the mean mapping quality is reported instead. Constants:
minimum clip 5 bases; clip position tolerance 5 bases; clip-fraction
window 50 bases; read-pair window 1000 bases; minimum span for read
pairs 1000 bases.

| Measure | Definition |
|---|---|
| `depth` | primary reads whose alignment covers POS |
| `mapq_mean`, `mapq0_frac` | mean mapping quality of those reads, and the fraction with mapping quality 0; empty when `depth` is 0 |
| `softclip_frac_w50` | among primary reads that overlap POS - 50 to POS + 50, the fraction with a soft clip of at least 5 bases at either end. The 3' clip of a read in a proper pair whose fragment (\|TLEN\|) is no longer than the read is not counted: that read runs into the adapter. Empty when no read overlaps |
| `clip_n` | reads (primary or supplementary) clipped, soft or hard, by at least 5 bases on the breakend's side, with the clip between POS + CIPOS lower - 5 and POS + CIPOS upper + 5. Orientation 1: the alignment ends there and its right end is clipped; -1: the alignment starts there and its left end is clipped |
| `split_n` | of the `clip_n` reads, those with an `SA` alignment on the mate's chromosome whose edge lies in the mate's window, by the same rule at the mate's POS, CIPOS and orientation (orientation 1: the SA alignment's last base; -1: its first base). Empty without a mate |
| `disc_n` | read pairs with a primary read on the breakend's side (orientation 1: forward strand, alignment start in POS - 1000 to POS; -1: reverse strand, alignment end in POS to POS + 1000) and its mate mapped on the mate breakend's side by the same rule (the mate's end from its `MC` tag, else its start plus the read length). Empty without a mate, and empty when both breakends are on one chromosome less than 1000 bases apart, where normal pairs fit the same geometry |
| `improper_frac_w1000` | among primary paired reads on the breakend's side (as for `disc_n`), the fraction whose pair the aligner did not flag proper (mate unmapped, on another chromosome, wrong orientation, or insert size outside the library's distribution). A measure of the local background, not of support. Empty when there is no such read |
| `support_mapq_mean` | mean mapping quality of the reads counted in `clip_n` and `disc_n`, each read once (its first record met); empty when there are none |

This is a recount from the alignments, not a copy of ESVEE's counts.
ESVEE assembles the reads, so it also counts reads that extend into an
inserted sequence or a repeat without a clip at POS; `clip_n` and
`split_n` can be lower than the VCF's `SF`. Measured on the NKTCL run
(2026-11-25): `bam_tumor_split_n` is 0 at both breakends for 4,940 of
the 10,145 two-breakend PASS events, all of which have `SF` of 1 or
more (750 have more than 50). A `split_n` of 0 in another BAM is
therefore not evidence that the event is absent there
(`SV_REPORT.md` section 6.7).

A task stops if the read groups of its BAM do not all name the sample
(`SM`), if a PURPLE VCF is missing, if a record lacks one ID, REF or
ALT, if an ID occurs twice in a VCF, if ALT neither starts nor ends
with REF, if a `MATEID` names no record, if a breakend's chromosome is
not in the BAM, or if the array has a different number of tasks than
there are BAMs.

**In the table.** For each breakend (`start_` or `end_` in front):

| Columns | Definition |
|---|---|
| `bam_tumor_<measure>` | the row of this record in the file of the group's tumor (`sample_id`) |
| `bam_normal_<measure>` | the same in the file of the group's normal (`normal_id`); empty for tumor-only groups |
| `bam_other_normals_<summary>`, `bam_other_tumors_<summary>` | over the normal BAMs, or the tumor BAMs, of all other groups, at this record: `n` the number of BAMs; `clip_samples`, `split_samples`, `disc_samples` the number of them with `clip_n`, `split_n`, `disc_n` above 0; `split_reads`, `disc_reads` the sums of `split_n` and `disc_n` |

In this cohort `bam_other_tumors_n` is 13 on every breakend, and
`bam_other_normals_n` is 8 for a matched group and 9 (every normal) for
a tumor-only group. The `split_*` summaries are empty where the
breakend has no mate, the `disc_*` summaries where `disc_n` is not
measured (the geometry is the same in every BAM).

`svCohort.R` stops if `--bam-dir` does not exist, if it holds a file
for a sample that is not in `samples.tsv`, if a tumor or normal has no
file, if a file's columns are not `group_id`, `id` and the nine
measures in that order, or if a file's rows do not match the cohort's
PURPLE records one to one. The joins to the events are one to one, and
`bam_tumor_depth` must be non-empty on every breakend.

htslib warns "The index file is older than the data file" for several
BAMs and one PURPLE SV VCF. In `out/` these index files have earlier
modification times than their data files (for example a 174 GB BAM:
index and the other files of its folder 13:40:06, BAM 13:56:25, all on
2026-11-18). The warning compares file times only; the indexes were
not rebuilt or compared with their files.

#### Population SV call sets

These columns say whether a breakend, or an event, matches a site of a
population SV call set. Four short-read WGS call sets are used; no
long-read set. The columns are counts and site IDs; they do not mark or
remove anything.

**Sources.** Downloaded by `fetchTracks.sh` on 2026-11-22 and listed
in `TRACK_DIR/tracks.tsv` with their full URLs. `md5 check` compares
the file with the md5 the source publishes, where it publishes one.

| Source | Call set | File | Records | md5 | md5 check |
|---|---|---|---|---|---|
| `gnomad_sv` | gnomAD-SV v4.1 sites | `gnomad.v4.1.sv.sites.vcf.gz` (Google Cloud bucket `gcp-public-data--gnomad`, `release/4.1/genome_sv/`) | 2,154,486 | 3ee614951c2f0c36659842876f7ce0ba | ok (the bucket's md5) |
| `kgp_sv` | 1000 Genomes high coverage, 3,202 samples (Byrska-Bishop 2022) | `1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz` (EBI FTP, `20210124.SV_Illumina_Integration/`) | 173,366 | 44be92aae63a812c19c2d6c91b22b4b9 | ok (the release manifest) |
| `dbvar_nstd166` | gnomAD-SV v2.1, 14,891 genomes (Collins 2020), as dbVar study nstd166 | `nstd166.GRCh38.variant_call.vcf.gz` (NCBI dbVar FTP) | 308,858 | 160623af53411e9b7b65d0826c980465 | none (dbVar lists no md5) |
| `dbvar_nstd200` | CCDG, 17,795 genomes (Abel 2020), as dbVar study nstd200 | `nstd200.GRCh38.variant_call.vcf.gz` (NCBI dbVar FTP) | 299,040 | 1f58bce51299daa768d81cb1469f0fcd | none (dbVar lists no md5) |

**Site tables, `scripts/popSvSites.sh`.** A job, once per site, after
`fetchTracks.sh`. `bcftools query` writes one line per VCF record,
every record kept (any FILTER), to `TRACK_DIR/<source>/sites.tsv.gz`:
`id`, `chrom`, `pos`, `end`, `chr2`, `end2`, `svtype`, `svlen`,
`filter`, `af` (VCF ID, CHROM, POS, INFO `END`, `CHR2`, `END2`,
`SVTYPE`, `SVLEN`, FILTER, INFO `AF`); a missing value is `.`. The
dbVar files have no `END2` and no FILTER (`.` on every line). For a
dbVar record whose ALT is in breakend notation (`N[chr:pos[` and the
like), `end2` is the partner position in ALT, and `chr2` its
chromosome where `CHR2` is missing. dbVar names the primary chromosomes
`1`-`22`, `X`, `Y`, `MT`; they are written `chr1`-`chr22`, `chrX`,
`chrY`, `chrM`, as in the run's GRCh38. Other contigs are kept as
written and match no breakend. The job stops if a table's line count
differs from the VCF's record count: from the index, or, for the 1000
Genomes file, whose index has no count metadata, from its data lines.

Site types and FILTER values in the tables (counted on 2026-11-23):

| Source | SVTYPE (records) | FILTER |
|---|---|---|
| `gnomad_sv` | DEL 1,206,278; BND 356,035; INS 304,645; DUP 269,326; CPX 15,189; INV 2,193; CNV 721; CTX 99 | PASS 1,199,117; other 955,369 |
| `kgp_sv` | DEL 90,259; INS 49,693; DUP 28,242; CPX 3,568; INV 920; CNV 673; CTX 11 | PASS on every record |
| `dbvar_nstd166` | DEL 154,952; INS 101,543; DUP 42,026; CNV 9,637; INV 700 | none |
| `dbvar_nstd200` | DEL 202,063; BND 51,125; DUP 44,149; INV 1,703 | none |

**Breakpoints of a site.** The first breakpoint is `chrom`:`pos`. The
second, where the site has one, is on `chr2` (or `chrom` where `chr2`
is missing) at `end2` (or `end` where `end2` is missing). In the
gnomAD table, DEL, DUP, INV, CNV and CPX records and all but 289 INS
records have no `END2`, so their second breakpoint is `end`; the CTX
records and 254,903 BND records have `CHR2` and `END2`. Sites carry no
orientation here, and none is compared; nor are confidence intervals
or inserted sequences.

gnomAD writes 101,132 of its BND records with `CHR2` = CHROM and no
`END2`. Their second breakpoint is `end`, which equals `pos` in 98,515
of them, so those sites have one breakpoint position. Prism does not
compute a partner position from `SVLEN`. All 101,132 have FILTER
`UNRESOLVED` (alone or with other values), so they can add to `_n`, not
to `_n_pass` or `_max_af`.

**PASS and AF.** A site is PASS when its FILTER is `PASS` or missing
(`.`), so every dbVar site is PASS. AF is the site's INFO `AF`.

**Per breakend.** For each breakend (the start of every event, the end
of two-breakend events) and each source:

| Column | Definition |
|---|---|
| `<source>_w100_n` | number of sites with a breakpoint on the breakend's chromosome and \|site position - POS\| <= 100; any type, any FILTER; each site counted once |
| `<source>_w100_max_af` | highest AF among those sites that are PASS; empty if there is none, or none has an AF |

**Per event, by breakpoints: `<source>_bp100_*`.** Two-breakend events
only. A site matches when one of its breakpoints is within 100 bases of
the event's start (same chromosome) and its other breakpoint within 100
bases of the event's end, either way round, and its SVTYPE is allowed
for the event's `info_SVTYPE`:

| Event `info_SVTYPE` | Site SVTYPE allowed |
|---|---|
| DEL | DEL, CNV, BND |
| DUP | DUP, CNV, BND |
| INV | INV, BND |
| INS | INS |
| BND | BND, CTX |

**Per event, by reciprocal overlap: `<source>_ro50_*`.** DEL, DUP and
INV events (both breakends on one chromosome). The event's interval
runs from the lower to the higher of its two positions. A site's
interval runs from the lower to the higher of `pos` and `end`; only
sites with an `end`, with `chr2` missing or equal to `chrom`, and of
an allowed type are used: for a DEL event DEL and CNV, for a DUP DUP
and CNV, for an INV INV. A site matches when the two intervals overlap
by at least half the length of the longer one (lengths in bases, both
ends included).

Each event rule gives four columns:

| Column | Definition |
|---|---|
| `_n` | number of matching sites, any FILTER |
| `_n_pass` | number of matching sites that are PASS |
| `_ids` | their site IDs (VCF ID), each once (a list) |
| `_max_af` | highest AF among the matching PASS sites |

`_n` and `_n_pass` are 0 where the rule applies and no site matches,
and empty where it does not apply. A site of another type (for example
CPX) never matches an event; it can still count in `<source>_w100_n`.

`svCohort.R` stops if a site table is missing, if its columns are not
the ten above in that order, or if an AF is neither one number nor
`.`.

#### Column dictionary, TSV format and manifest

`sv_cohort_filtered.columns.tsv` has one row per column of the table,
in table order:

| Column | Content |
|---|---|
| `column` | column name |
| `description` | source and meaning; INFO and FORMAT descriptions from the VCF header (the first group's text for each key); notes for split parts, PURPLE's per-breakend values and event-level columns; the `xs_*` descriptions end with the definition of the sets |
| `non_empty` | number of non-empty cells |
| `type` | `string`, `integer`, `double`, `boolean` or `list<...>` (see Types) |
| `list_sep` | separator of the list in the TSV (`\|`, `,` or `;`); empty for scalars |
| `list_count_column` | the column with the number of elements of each list, where one exists |

The run stops if a column has no description.

A reference copy of the dictionary is kept in the repo,
`docs/sv_cohort_filtered.columns.csv` (CSV, for people to read),
without `non_empty`, which depends on the cohort. Everything else in
it is fixed by the code, the ESVEE and PURPLE versions (INFO and FORMAT
descriptions and types come from the VCF headers) and `tracks.tsv`.
A different cohort can still add or remove columns: an `info_*`
column exists only for an INFO key that occurs in the data, and
`normal_*` columns only if a group has a matched normal.
`postOA.sh` passes the copy as `--dict-ref`. Before the table is
written, `svCohort.R` compares the new dictionary with it. Any added,
removed or reordered column, or any changed description, type,
`list_sep` or `list_count_column`, stops the run with the list of
differences, and the run writes no table. The new dictionary is then
written as `sv_cohort_filtered.columns.new.csv` in the database
folder. If the change is intended, copy that file over the reference,
commit it and rerun, so every column change shows up as a reviewed
commit. The comparison was tested on the current build (no
differences) and on one edit of each kind: a changed description,
type and `list_sep`, an added and a removed column, and a reordering.
Each edit was reported.

The TSV is written without quoting, with an empty cell for a missing
value, so every reader sees the values as they are. A value with a tab,
a newline or a double quote would need quoting; if one is found, the
run stops.

The table (TSV and Parquet) and the dictionary replace their rows in
`manifest.tsv`.

#### When the script stops

`svCohort.R` stops, and writes no table, when:

1. an argument is missing or unknown;
2. the PURPLE SV VCF has no column named `sample_id`, or more than two
   sample columns;
3. a VCF ID occurs twice in a group;
4. a record other than FILTER `INFERRED` has no ESVEE row with the same
   ID, chromosome, position and orientation;
5. `MATEID` links are not reciprocal;
6. the events do not account for every VCF record (2 x two-breakend +
   one-breakend events = records);
7. an INFO or FORMAT key has a different `Number` or `Type` in the VCF
   headers of two groups;
8. a LINX `vcfIdEnd` is neither the event's end ID nor its start ID;
9. a LINX `svs` row matches no event;
10. LINX's start is the VCF end for an SV whose breakends differ in
    position;
11. a LINX `breakend.tsv` row has an `svId` of no LINX event or an
    `isStart` other than `true`/`false`, or the rows placed on event
    ends do not number exactly the rows of `breakend.tsv`;
12. a value to be joined with `|` already contains `|`;
13. a row breaks the assumption of an unpacking rule (rules 1 to 5);
14. the number of columns changes when they are put in order;
15. a column has no description, or a VCF column has no header
    definition;
16. a value contains a tab, a newline or a double quote;
17. the column dictionary differs from its reference copy
    `docs/sv_cohort_filtered.columns.csv` (see above);
18. the driver gene panel is missing, lists a gene twice or has a
    `likelihoodType` other than ONCO or TSG, or a `gene_role` list does
    not line up with `linx_bnd_gene`;
19. a fusion breakend ID has no `breakend.tsv` row, or a fusion name is
    neither `<5' gene>_<3' gene>` nor `@<region>_<3' gene>`;
20. `--bam-dir` does not exist, holds a file for a sample not in
    `samples.tsv`, or lacks the file of a tumor or normal; a file's
    columns are not the expected ones; or a file's rows do not match
    the cohort's PURPLE records one to one (see BAM evidence);
21. a population site table is missing, has other columns, or has an
    AF that is neither one number nor `.` (see Population SV call
    sets).

Not checked by the script: that `drivers.tsv` cluster IDs are in `clusters.tsv`
(`docs/scripts/checkLinx.R` checks it); that `MATEID` is reciprocal
in the cross-sample sets; that `sample_id` is unique over groups (which
`event_key` assumes); that a population site ID contains no `|`.

The computed columns (region tracks, cross-sample counts, BAM evidence,
population matches) are not recomputed by step 4, which rebuilds the
source values. Section 1.5 recomputes them for a random sample of
events by a second method.

### 1.3 Step 3, `scripts/tsvToParquet.py`: the typed Parquet copy

```
python3 Prism/scripts/tsvToParquet.py TABLE.tsv.gz DICTIONARY.tsv OUT.parquet
```

`postOA.sh` runs it on the table and dictionary of step 2 and writes
`sv_cohort_filtered.parquet` (pyarrow, zstd compression). It reads the
TSV with every column as string (no quote handling, an empty cell is
null) and converts each column to the type its dictionary row gives:

- `integer`: 32-bit integers, so R reads them as integer (R has no
  64-bit integer type); a value outside that range stops the run.
- `double`: 64-bit floating point (`NaN` and `Inf` are accepted).
- `boolean`: `true`/`TRUE` is true, `false`/`FALSE` is false.
- `list<T>`: split on `list_sep`, each element converted to T, an empty
  element null. Where the dictionary names a `list_count_column`, a row
  with a count above 0 always has a list (so a list of one empty element
  is kept), and every list must have exactly that many elements.

Every value is first matched against the pattern of its type; any value
that does not fit stops the run with the column and an example value.
The table and dictionary columns must be the same, in the same order.
After writing, it reads the file back and checks the number of rows,
the column names and, for every scalar column, the number of nulls. It
exits with status 1 on any error.

**Reading the table.** Use the Parquet file; it carries the types and
the lists. In R, start from `Prism/examples/readSvCohort.R`: source it
from the project folder (or set `SV_DB_DIR` first) and it loads the
table as `sv` and the dictionary as `sv_dict`, checks that they have
the same columns, and defines `sv_find()` (dictionary rows whose name
or description matches a pattern) and `sv_transcripts()`,
`sv_fusions()`, `sv_links()`, `sv_cluster_drivers()` (one row per
element of a list family, with the event's keys in front), and
`sv_flatten()` (one row per event, every list column joined into one
string with `|`; a missing element is an empty slot, so positions stay
aligned; an empty list is NA). Its comments
end with worked examples.

```r
source("Prism/examples/readSvCohort.R")
pass <- sv |> filter(filter == "PASS")
sv_transcripts(pass) |> filter(canonical, disruptive) |> distinct(group_id, gene)
```

Without the bootstrap, in R (nanoparquet; list columns become R lists):

```r
library(tidyverse)
sv <- nanoparquet::read_parquet("res/prism/db/sv_cohort_filtered.parquet")
dict <- read_tsv("res/prism/db/sv_cohort_filtered.columns.tsv", show_col_types = FALSE)

# Genes at the start breakend of PASS events in LINX, canonical transcripts:
# unnest the aligned lists of one family together.
sv |>
  filter(filter == "PASS", in_linx) |>
  select(event_key, start_linx_bnd_gene, start_linx_bnd_canonical) |>
  unnest(c(start_linx_bnd_gene, start_linx_bnd_canonical)) |>
  filter(start_linx_bnd_canonical)
```

`unnest()` drops rows whose list is empty; add `keep_empty = TRUE` to
keep them. In Python (pyarrow, with pandas or polars):

```python
import pyarrow.parquet as pq
sv = pq.read_table("res/prism/db/sv_cohort_filtered.parquet").to_pandas()
genes = sv[["event_key", "start_linx_bnd_gene", "start_linx_bnd_canonical"]].explode(
    ["start_linx_bnd_gene", "start_linx_bnd_canonical"])
```

The TSV holds the same values as text, for tools without Parquet
support. Read it with every column as character, no quote handling, and
an empty cell as missing (in R:
`read_tsv(f, col_types = cols(.default = "c"), quote = "", na = "")`),
then use the dictionary's `type` and `list_sep` to convert.

### 1.4 Step 4, `docs/scripts/checkSvCohort.R`: validation

```
Rscript --vanilla Prism/docs/scripts/checkSvCohort.R [OUT_DIR] [DB_DIR] [HMF_DIR]
```

`HMF_DIR` (for the driver gene panel) defaults to the Hartwig resource
folder of the run, as in `postOA.sh`; `postOA.sh` passes it.

`postOA.sh` runs it as its last step and fails if any check fails (exit
status 1). It can also be run on its own after a change to `svCohort.R`
or a pipeline upgrade. It undoes every unpacking rule: it joins split
pairs, rebuilds PURPLE's two-value fields from both breakends, copies
event-level values to both records, and rebuilds each dropped column
from the column it repeats. From that it rebuilds every VCF record and
every ESVEE row, and compares them with the source files (a VCF `.`
compares equal to empty).

Per group, each check must be TRUE:

| Check | What it compares |
|---|---|
| `records_once` | every record of `purple.sv.vcf.gz` is the start or the end of exactly one event, and no other ID is in the table |
| `fields_equal` | CHROM, POS, REF, ALT, QUAL, FILTER of each rebuilt record equal the VCF |
| `info_equal` | the rebuilt INFO key=value pairs equal the VCF INFO of each record (flags as `true`) |
| `format_equal` | the rebuilt FORMAT values of the tumor and the normal equal the VCF |
| `esvee_equal` | every column of every ESVEE `breakend.tsv` row of a record equals the rebuilt value (StrandBias: FORMAT SB of each sample, within 0.005); `INFERRED` records have no ESVEE row |
| `linx_ok` | events in LINX = LINX `svs.tsv` rows, and every event with `linx_start_is_vcf_start` FALSE has both breakends at one position |
| `linx_bnd_coords_ok` | every LINX breakend row, placed by `isStart`, has the coordinates of the end it is placed on (end rows of one-breakend events are counted, not compared: LINX writes them at a location that is not a VCF record) |

Over the whole table:

| Check | What it compares |
|---|---|
| `linx_tables_equal` | every LINX value in the table, as text, against the stacked LINX tables (`data/linx_somatic_*.tsv.gz`), each source row exactly once and nothing extra: the `linx_*` columns against `svs`; the `linx_cluster_*` columns, once per cluster, against `clusters`; the per-transcript rows rebuilt from the start and end lists against `breakend`; the fusions from their 5' side and from their 3' side, each on the SV of LINX's breakend ID, against `fusion`; the links from the lower and from the upper SV against `links`; the drivers of each event's cluster against `drivers` (the same list on every event of a cluster; drivers with no SV cluster are counted, not in the table) |
| `gene_roles_equal` | every element of `start_gene_role` and `end_gene_role` equals the panel's `likelihoodType` of the gene in the same position of `linx_bnd_gene` (empty if not on the panel), and the list lengths equal `linx_bnd_n` |
| `fusion_genes_equal` | every element of `fusion_gene_5p` and `fusion_gene_3p` equals the `breakend.tsv` gene of the fusion's breakend ID in the same position, and every `linx_fusion_name` is `<5' gene>_<3' gene>` or `@<region>_<3' gene>` |
| `parquet_equal` | every cell of the Parquet file equals the TSV cell converted to the dictionary type, by a separate implementation in R (nanoparquet) of the conversion in `tsvToParquet.py`; lists compared element by element. Doubles are compared to 1 part in a million (`same_values()`): R's `as.numeric` does not always give the nearest double (`5.754e-03` lands one unit in the last place from the value pyarrow and Python give, checked with exact arithmetic). Conversion errors, which this check is for, are far larger; every other type must be identical |

It also prints the dictionary's type counts, the columns empty in every
row, the events by number of breakends, LINX and FILTER, and counts of
the region and cross-sample columns.

**Results.**

- 2-group test, 2026-11-22: two groups of this project (one with a
  matched normal, one tumor-only; a reduced track folder of 18 tracks
  and part of RepeatMasker), 2,025 VCF records, 1,109 events, 341
  columns. Every check passed, including `parquet_equal`. Column types:
  string 48, integer 122, double 19, boolean 47, `list<string>` 55,
  `list<integer>` 34, `list<double>` 4, `list<boolean>` 12. One column
  was empty in every row (`end_info_INFERRED`, as inferred breakends are
  always one-breakend events).
- Full cohort, 2026-11-22 (`postOA.sh` job 13234922, all 20 tracks;
  log `out/logs/other/postOA.13234922.out`): 14 groups (9 with a
  matched normal, 5 tumor-only), 32,965 VCF records, 21,890 events, 346
  columns. 61 duplicate columns were dropped and 36 columns are
  event-level. Every per-group check was TRUE in all 14 groups, and
  `parquet_equal` was TRUE. Column types: string 48, integer 122,
  double 19, boolean 52, `list<string>` 56, `list<integer>` 33,
  `list<double>` 4, `list<boolean>` 12. One column was empty in every
  row (`end_info_INFERRED`). `linx_tables_equal` was added afterwards
  and run on this table (job 13239483, log
  `out/logs/other/checkSvCohort.13239483.out`): every LINX table check
  TRUE; 23 LINX driver rows belong to no event's cluster and are not in
  the table. A deliberately changed gene name in a copy of the 2-group
  table made the check fail on that row. The cohort includes the two cases the
  2-group test lacks: SVs with both breakends at identical coordinates
  (4 INV; their LINX breakend rows are each placed once, on the end
  given by `isStart`) and 4 single-breakend events with 5 LINX end
  rows. 65 events had LINX's start at the VCF end, all with both
  breakends at one position, all in tumor-only groups. Events by number
  of breakends, LINX and FILTER:

  | Breakends | In LINX | FILTER | Events |
  |---|---|---|---|
  | 1 | no | PON | 503 |
  | 1 | yes | INFERRED | 8,241 |
  | 1 | yes | PASS | 2,071 |
  | 2 | no | PON | 930 |
  | 2 | yes | PASS | 10,145 |

- Full cohort with BAM evidence and population SV call sets,
  2026-11-23 (`postOA.sh` job 13303561, the code of commit 7786c4d;
  log `out/logs/other/postOA.13303561.out`): the same 21,890 events,
  454 columns. The 108 new columns: 60 BAM columns (30 per breakend),
  16 per-breakend population columns (8 per breakend) and 32 per-event
  population columns. 23 BAM files (14 tumors, 9 normals) of 32,965 rows each,
  758,195 rows. The dictionary equaled its reference copy. Every check
  was TRUE, including `linx_tables_equal` and `parquet_equal`. Column
  types: string 48, integer 186, double 55, boolean 52, `list<string>`
  64, `list<integer>` 33, `list<double>` 4, `list<boolean>` 12. One
  column was empty in every row (`end_info_INFERRED`). The BAM
  columns: `bam_tumor_split_n` is non-empty on the 11,075 two-breakend
  events, `bam_tumor_disc_n` on 7,821 of them (the other 3,254 have
  both breakends on one chromosome less than 1000 bases apart), and
  `start_bam_normal_*` on the 15,739 events of matched groups.
  Population matches (of 32,965 breakends, 11,075 two-breakend events
  and 4,685 DEL, DUP and INV events):

  | Source | Sites | Breakends with a site within 100 bases | Events matched by breakpoints | Events matched by reciprocal overlap |
  |---|---|---|---|---|
  | `gnomad_sv` | 2,154,486 | 15,392 | 1,307 | 1,347 |
  | `kgp_sv` | 173,366 | 3,283 | 478 | 559 |
  | `dbvar_nstd166` | 308,858 | 3,494 | 444 | 556 |
  | `dbvar_nstd200` | 299,040 | 2,912 | 492 | 635 |

- Full cohort with the lead block, gene roles and fusion genes,
  2026-11-23 (`postOA.sh` job 13347029, branch `feat/sv-columns`; log
  `out/logs/other/postOA.13347029.out`): the same 21,890 events, 458
  columns (454 + `start_gene_role`, `end_gene_role`, `fusion_gene_5p`,
  `fusion_gene_3p`), the first 53 in the order of the lead block. The
  first run (job 13346246) stopped at the dictionary check as intended
  and listed the 4 added columns; the new dictionary was reviewed (458
  columns; no description, type or list field of an existing column
  changed; existing columns only reordered) and copied over
  `docs/sv_cohort_filtered.columns.csv`. Every check was TRUE,
  including `gene_roles_equal`, `fusion_genes_equal`,
  `linx_tables_equal` and `parquet_equal`. Column types as in the
  previous build, with `list<string>` 68. The BAM, population, region
  and cross-sample counts equal those of job 13303561.
  - Gene role elements: start breakends ONCO 345, TSG 752, empty
    42,485; end breakends ONCO 198, TSG 449, empty 20,767.
  - Fusion elements: 2,342 (1,171 fusions, each listed once for the SV
    of its 5' breakend and once for the SV of its 3' breakend: 1,149
    fusions have both on one SV, which lists them twice; 22 chained
    fusions are on two SVs, once on each); names `<5' gene>_<3' gene>`
    2,318, enhancer-named 24 (the 12 fusions in Joins, each made by one
    SV and listed twice on it), other 0. 9 distinct reported fusions
    (group, name).

  These counts are the ones `SV_Cohort_Code_Changes.md` (project
  folder) predicted from the previous table. The spot check (section
  1.5) was run on this table too (job 13361310), with every value
  equal.

- Rebuild from scratch, 2026-11-23, after `res/` was removed: `postOA.R`
  alone (job 13376213) to write `samples.tsv`, `bamEvidence.sh` (array
  job 13376222, 23 tasks, 23 BAMs, 758,195 rows), then `postOA.sh` (job
  13376226, log `out/logs/other/postOA.13376226.out`, Prism `180d01f`
  with no uncommitted changes). Result: 21,890 events x 458 columns,
  column dictionary equal to the reviewed copy, every check TRUE
  (`linx_tables_equal`, `gene_roles_equal`, `fusion_genes_equal`,
  `parquet_equal`, all per-group checks); one column empty in every row
  (`end_info_INFERRED`), as before. This is the table the SV report
  (section 2) reads.

#### Changes after review

Found while writing this section and fixed in `svCohort.R` on
2026-11-22:

- The dictionary described `region_hmf_sv_prep_blacklist` as "regions
  ESVEE prep leaves out"; no task of this run reads that file (above).
  The description now says so.
- `sv_prep_blacklist.38.bed` was read as 1-based; it is now read as
  Hartwig's loader reads it, [PosStart + 1, PosEnd]. Only a breakend at
  exactly `PosStart` changes.
- The dictionary text for `region_centromere_distance` said "0:
  inside"; it now also says that 0 means the base next to an interval,
  and that the cell is empty with no model on the chromosome.
- `own_germline_*` was empty, not 0, for a matched group whose
  germline file has no records; matched groups are now taken from
  `normal_id` in `samples.tsv`.

### 1.5 `docs/scripts/spotCheckSvCohort.R`: spot check of the computed columns

```
Rscript --vanilla Prism/docs/scripts/spotCheckSvCohort.R OUT_DIR DB_DIR TRACK_DIR HMF_DIR [N] [N_BAM]
```

Run it by hand after a build, as a job (about 10 minutes for the
defaults); `postOA.sh` does not run it. It needs `bcftools`
(`BCFTOOLS`, default `/admin/software/bcftools/bcftools-1.20/bin/bcftools`)
and `samtools` on the `PATH`.

Step 4 rebuilds the source values from the table; it does not
recompute the columns Prism computes. This script recomputes them for
a random sample of events by a second method, written separately from
`svCohort.R` (a plain filter per breakend instead of interval
overlaps; `bcftools` and `samtools` instead of the R and pysam
readers), and compares the result with the TSV. It samples N events
(default 300, any FILTER, fixed seed 20261122), takes their breakends
(the start, and the end of two-breakend events), and of those N_BAM
breakends (default 40) for the BAM recount. It exits with status 1 if
any value differs, and prints the differences.

| Check | Columns | Second method |
|---|---|---|
| `region_bed` | the 17 GIAB tracks and the ENCODE blacklist (18) | the BED rows with start < POS <= end |
| `region_labels` | ENCODE blacklist reason; RepeatMasker name, class, family (4) | the set of labels of the intervals that contain POS |
| `region_hmf` | `region_hmf_sv_prep_blacklist`, `region_hmf_unmap` (2) | the two files, with 1 added to the start of the `.bed` file only |
| `centromere` | `region_centromere_distance` (1) | the bases strictly between POS and the nearest centromere model interval of the chromosome (0 inside or adjacent) |
| `xs_breakend` | `xs_<set>_groups_exact`, `_w50`, `_group_list_w50`, `own_germline_*` (11) | every record of the three sets, read with `bcftools query` |
| `xs_pair` | `xs_pair_<set>_*` (6) | every `MATEID`-linked pair of the sets |
| `bam_own` | `bam_tumor_depth`, `bam_tumor_clip_n`, `bam_tumor_split_n` (3) | `samtools view -F 0x704` of the breakend's region, the SAM text parsed in R |
| `bam_others` | `bam_other_normals_*`, `bam_other_tumors_*` (12) | recomputed from the per-BAM files |
| `pop_breakend` | `<source>_w100_n`, `_max_af` (8) | every site of the source table |
| `pop_event` | `<source>_bp100_*`, `<source>_ro50_*` (32) | every site of the source table |

What it does not cover: of the BAM measures, only `depth`, `clip_n`
and `split_n` of the own tumor are recounted, on N_BAM breakends;
`bam_others` and the population checks start from the same per-BAM
files and site tables as `svCohort.R`, so they check the summaries and
the matching, not the recount (`bamEvidence.py`) or the reduction to
sites (`popSvSites.sh`, checked by its record counts). Lists are
compared as sets.

**Result.** Job 13317208, 2026-11-23, on the table of job 13303561
(log `spotCheck.13317208.out` in the project folder): 300 events, 458
breakends, 40 of them recounted with samtools. Every value was equal:
31,772 values in 97 columns.

Job 13361310, 2026-11-23, on the 458-column table of job 13347029
(log `spotCheck.13361310.out`): the same seed, 300 events and 458
breakends, and the same result, every value equal (31,772 values in 97
columns; per-check counts as below). The columns added in job 13347029
(gene roles, fusion genes) are checked by `checkSvCohort.R`, not here.

| Check | Columns | Values | Equal |
|---|---|---|---|
| `region_bed` | 18 | 8,244 | 8,244 |
| `region_labels` | 4 | 1,832 | 1,832 |
| `region_hmf` | 2 | 916 | 916 |
| `centromere` | 1 | 458 | 458 |
| `xs_breakend` | 11 | 5,038 | 5,038 |
| `xs_pair` | 6 | 948 | 948 |
| `bam_others` | 12 | 5,496 | 5,496 |
| `bam_own` | 3 | 120 | 120 |
| `pop_breakend` | 8 | 3,664 | 3,664 |
| `pop_event` | 32 | 5,056 | 5,056 |

Two earlier runs failed because of errors in the spot-check script, not
in the table. Job 13304919 stopped when a helper was given integers
(it now converts to text first). Job 13314940 reported 11 of 40
`bam_tumor_split_n` values as 0 where the table is empty: the script
read the empty `end_id` of one-breakend events as a mate ID. Now it
reads it as missing, and the table's value (empty without a mate) is
confirmed.

## 2. `scripts/svReport.R`: the SV report

The SV report is the workbook people read about the somatic SVs:
one row per PASS event with the genes and fusions LINX gives it, five
evidence flags, and a high-confidence set (events with no flag), plus
gene and fusion sheets that count tumors over high-confidence and over
all PASS events. It reads only the SV table of section 1
(`res/prism/db/sv_cohort_filtered.parquet`) and writes
`res/prism/sv/prism_sv_report.xlsx` and `prism_sv_events.csv`:

```
Rscript --vanilla Prism/scripts/svReport.R [DB_DIR] [OUT_DIR]
```

Its full description is a document of its own, `SV_REPORT.md`: why it
exists, every sheet and column, the five flags with their exact
definitions and the SV table columns behind them, the LINX rules behind
the gene and fusion columns (checked in the LINX code at `linx-v2.3.1`),
the checks, the issues met building it, and the results of this run.
Keep the two in step: a change to `svReport.R` changes `SV_REPORT.md`.

Two exceptions to the project rule, both user decisions of 2026-11-23
and both limited to summary sheets (`events_all` keeps every PASS event
and every LINX fusion of it):

- LINX fusions whose 5' and 3' breakends are in the same gene join two
  parts of one gene (intragenic rearrangements, which the user treats
  as cryptic splice products), not two genes. They are left out of the
  fusion sheet, named `fusions_real`, and of the gene fusion counts
  (96 fusions on 101 events in this run, counted in the `rule` sheet).
- The `genes` sheet lists genes with 2 or more tumors in at least one
  count.

The high-confidence set follows the project rule: it is a column
(`high_confidence`, with `n_flags` and the five `flag_*` columns), and
the `events_all` sheet keeps every PASS event. The flag rule is the
user's decision of 2026-11-23.

**Result.** 2026-11-23, on the table of `postOA.sh` job 13376226:
12,216 PASS events, 4,180 high confidence (no flag); flags set:
LINX ARTIFACT 4,233, region 3,363, other patients 2,338, population
1,006, own normal 34 (`SV_REPORT.md` section 6.6).

**Review.** 2026-11-25: every sheet recomputed from the TSV copy of the
table with separate code and compared with the workbook, no cell
different. Checks added to `svReport.R` for a missing column, a missing
FILTER value, the LINX transcript lists and the fusion entries
(`SV_REPORT.md` section 9); the workbook is unchanged. The same review
measured what "no flag" does not cover (`SV_REPORT.md` section 6.7):
the BAM recount finds no split read in the event's own tumor for 4,940
of the 10,145 two-breakend PASS events; 1,051 high-confidence
two-breakend events share one breakend with another group's calls; and
1,395 of the 1,883 high-confidence translocations are in the four
groups PURPLE fits at purity 0.18 to 0.21.

## 3. `scripts/orangeSvReport.R`: the ORANGE SV report

ORANGE's JSON (`out/<group>/orange/<sample>.orange.json`) lists the
pipeline's own short list of SVs: LINX fusions with `reported` true,
and somatic and germline LINX breakend rows with `reportedStatus`
REPORTED (ORANGE 5.0.1 source, `LinxDataLoader`,
`LinxBreakendInterpreter`). ORANGE tests no SV FILTER value; its SVs
are the ones LINX analyzed (somatic PASS and INFERRED, germline PASS
and PON). This script collects those records of every sample, with the
ORANGE fields unchanged, places each on its SV table event (by sample
and LINX `svId`; fusions through their LINX `fusion.tsv` row) and adds
the event's FILTER and its SV report flags and `high_confidence`. It
writes `res/prism/sv/orange_sv_report.xlsx`:

```
Rscript --vanilla Prism/scripts/orangeSvReport.R [OUT_DIR] [DB_DIR] [SV_DIR]
```

Run it after `svReport.R`. It stops if the JSON records are not exactly
the LINX rows the ORANGE code selects, if a record cannot be placed on
its event, or if `prism_sv_events.csv` does not hold the PASS events of
the SV table. Full description: `SV_REPORT.md` section 13.

**Result.** 2026-11-25, on the table of `postOA.sh` job 13376226: 9
fusions on 9 PASS events (3 high confidence), 159 somatic breakend
records on 152 PASS events (71 high confidence), 2 germline breakends;
no record on an INFERRED or PON event.
