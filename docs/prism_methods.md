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

## Software

| Software | Version | Source of the version |
|---|---|---|
| nf-core/oncoanalyser | 3.0.0 (`v3.0.0-g7c74c87`), Nextflow 25.10.4 | `out/pipeline_info/nf_core_oncoanalyser_software_mqc_versions.yml` |
| ESVEE | 2.0.1 | same file |
| PURPLE | 4.4 | same file |
| LINX | 2.3.1 | same file; also `samples.tsv` column `linx_version` |
| R | 4.5.1 | `Rscript` on the `PATH` of the job |
| Bioconductor | 3.22: GenomicRanges 1.62.1, IRanges 2.44.0, S4Vectors 0.48.1 | same |
| CRAN | tidyverse 2.0.0 (dplyr 1.2.1, readr 2.2.0, tidyr 1.3.2, stringr 1.6.0, purrr 1.2.2), fs 2.1.0, nanoparquet 0.5.2 | same |
| Python | 3.10.14, pyarrow 25.0.1 | `python3` on the `PATH` of the job |

The versions above were read on 2026-11-22 from the `Rscript` and
`python3` that `postOA.sh` finds on the `PATH`. `postOA.sh` does not log
them. To record them for a run, use the same environment as the job and
run (a few seconds):

```
Rscript --vanilla -e "cat(R.version.string, '\n'); for (p in c('BiocVersion', 'GenomicRanges', 'IRanges', 'S4Vectors', 'tidyverse', 'dplyr', 'readr', 'tidyr', 'stringr', 'purrr', 'fs', 'nanoparquet')) cat(p, as.character(packageVersion(p)), '\n')"
python3 -c "import sys, pyarrow; print(sys.version.split()[0], 'pyarrow', pyarrow.__version__)"
```

The Prism version of a run is the `PRISM=` line at the top of its log
(`git describe --always --dirty` of the Prism folder).

## 1. postOA.sh: cohort tables

```
./Prism/scripts/fetchTracks.sh [TRACK_DIR]           # once per site, login node
sbatch ./Prism/postOA.sh [-t TRACK_DIR] [OUT_DIR]
```

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

Every file is read as text (`col_character`, no quote handling, no
value turned into NA), so values reach the table exactly as written.
The one exception is the VCF missing value `.`: in QUAL, FILTER, INFO
and FORMAT it is stored as an empty cell.

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

1. key columns (table below);
2. event-level columns: fields with the same value on both breakends,
   stored once, without a `start_`/`end_` prefix (see Unpacking, rule
   5);
3. every `start_*` column, then every `end_*` column (the same set of
   columns for each breakend);
4. the per-event LINX columns (`linx_*`);
5. the `xs_pair_*` columns.

Within a breakend, `info_*`, `tumor_*`, `esvee_bnd_*` and `normal_*`
columns appear in the order their keys are first seen, group by group,
and a split pair takes the place of its field. So the order can differ
between runs. Select columns by name.

Key columns:

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
| `region_*` | region tracks (see Region tracks) |
| `xs_*`, `own_germline_*` | the same breakend in other files (see Cross-sample columns) |

Per event: `linx_<column>` (LINX `svs.tsv`), `linx_cluster_<column>`
(`clusters.tsv`), `linx_fusion_*`, `linx_link_*`,
`linx_cluster_driver_*` (lists), and `xs_pair_*`.

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
  (`5p`, `3p`).
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

**Empty cells.** In the TSV an empty cell is a missing value; in the
Parquet copy it is null. A cell is empty when: the column is an
`end_*` column of a one-breakend event (except LINX end rows, see
Joins); the record has no such INFO or FORMAT key, or the VCF has `.`;
the group is tumor-only (`normal_*`, `own_germline_*`); the record is
`INFERRED` (`esvee_bnd_*`); the event is not in LINX (`linx_*`); there
is no source row for a list; or there is no interval or group to list.

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
   and `own_germline_*`: integer; group lists `list<string>`.
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

The 20 tracks, as downloaded by `scripts/fetchTracks.sh` on 2026-11-21
(`TRACK_DIR/tracks.tsv` holds the same list with full URLs). `md5 check`
is the comparison with GIAB's published md5 list
(`GRCh38-genome-stratifications-md5s.txt`); the other sources publish
none, so their md5 identifies the file used.

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
16. a value contains a tab, a newline or a double quote.

Not checked by the script: that fusion breakend IDs are in
`breakend.tsv` and that `drivers.tsv` cluster IDs are in `clusters.tsv`
(`docs/scripts/checkLinx.R` checks both); that `MATEID` is reciprocal
in the cross-sample sets; that `sample_id` is unique over groups (which
`event_key` assumes).

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
element of a list family, with the event's keys in front). Its comments
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
Rscript --vanilla Prism/docs/scripts/checkSvCohort.R [OUT_DIR] [DB_DIR]
```

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
