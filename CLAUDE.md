# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Prism is a wrapper for running nf-core/oncoanalyser 3.0.0 (mode `wgts`,
genome `GRCh38_hmf`) on IRIS (SLURM, Singularity). It holds no pipeline
code. Nextflow fetches the pipeline; Prism supplies the run script,
the configs and the samplesheet builder.

Prism is cloned into a project folder and every command is run from
that project folder, not from inside `Prism/`. Outputs (`out/`,
`.nextflow/`, `.nextflow.log`, `oncoanalyser.<jobid>.out`) go to the
project folder. Nothing a run writes belongs in the repo.

Read `docs/NF_SCRIPTING_INFO.md` before changing the run script or any
config. It lists the rules these files follow and the mistakes that led
to each rule. It is also the starting point for the next pipeline
wrapper, which gets its own repo and a copy of that file.

## Commands

Run from the project folder (the parent of `Prism/`):

```
./Prism/00.SETUP.sh                                     # once per clone: install Nextflow in Prism/bin
Rscript Prism/scripts/makeOncoanalyserInput.R mapping.csv manifest.csv oncoanalyser_input.csv
sbatch ./Prism/runOncoanalyser.sh INPUT.csv [WORK_DIR]  # as a SLURM job
./Prism/runOncoanalyser.sh INPUT.csv [WORK_DIR]         # in the current shell
./Prism/deliverLocal.sh [OUT_DIR] [RESULTS_DIR]         # hard link reports into results/{PDF,HTML,PNG}
sbatch ./Prism/qcCohort.sh [OUT_DIR]                    # QSEE over all groups into res/prism/qc/hmf/
sbatch ./Prism/cohortReport.sh [-m MANIFEST] [-g GENES] [OUT_DIR]  # cohort report into res/prism/cohort/
./Prism/scripts/fetchTracks.sh [TRACK_DIR]              # once per site: region tracks for postOA.sh (login node, downloads only)
sbatch ./Prism/scripts/popSvSites.sh [TRACK_DIR]        # once per site, after fetchTracks.sh: reduce the population SV VCFs to sites.tsv.gz
sbatch --array=0-N ./Prism/bamEvidence.sh [OUT_DIR]     # per run, before postOA.sh: read counts at every SV breakend in every BAM (./Prism/bamEvidence.sh prints N)
sbatch ./Prism/postOA.sh [-t TRACK_DIR] [OUT_DIR]       # cohort tables into res/prism/db/: stacked LINX tables, SV event table (Parquet, TSV), checked
Rscript --vanilla Prism/docs/scripts/checkLinx.R [OUT_DIR]  # recheck the LINX rows of docs/ACTUAL_METHODS.md
Rscript --vanilla Prism/docs/scripts/checkSvCohort.R [OUT_DIR] [DB_DIR] [HMF_DIR]  # rebuild the sources from sv_cohort_filtered and compare; recompute gene roles and fusion genes; check the Parquet copy (postOA.sh step 4)
Rscript --vanilla Prism/docs/scripts/spotCheckSvCohort.R OUT_DIR DB_DIR TRACK_DIR HMF_DIR  # recompute the derived columns of a random sample by a second method (as a job: about 10 minutes)
Rscript --vanilla Prism/scripts/svReport.R [DB_DIR] [OUT_DIR]  # SV report workbook (PASS events, LINX genes and fusions, five flags, high-confidence set) into res/prism/sv/; about a minute
Rscript --vanilla Prism/scripts/orangeSvReport.R [OUT_DIR] [DB_DIR] [SV_DIR]  # after svReport.R: ORANGE JSON SV records (fusions, somatic and germline breakends) of all samples, linked to SV table events and SV report flags, into res/prism/sv/orange_sv_report.xlsx; seconds
Rscript Prism/docs/dbScheme/scripts/fetchOrangeDatamodel.R orange-v5.0.1 DATAMODEL_DIR  # ORANGE data model sources
Rscript Prism/docs/dbScheme/scripts/buildScheme.R --datamodel DATAMODEL_DIR --out-dir out  # rebuild docs/dbScheme tables
```

`docs/ACTUAL_METHODS.md` records what each step actually does, checked
against code, `out/logs/<group>/<task>.command.{sh,log}` and output
files; tool READMEs are not trusted. Keep its status labels
(CONFIRMED, README, VARIANCE, OPEN) and recheck after any upgrade.

`docs/dbScheme/` records, field by field, which output file first holds
each per-sample value (role primary) and which files only copy or
derive it (ORANGE, QSEE, CUPPA `pred_summ`, Virus Interpreter), as the
basis for the cohort aggregation (`postOA.sh`).
`curated/` is edited by hand; `tables/` and `dbScheme.xlsx` are
generated from it, the ORANGE data model and a finished `out/`. What
ORANGE computes is in its source code, not its data model; read
`orange/src/main/java/` again on an ORANGE upgrade. Prism is public:
the generated tables hold only all/some/none and counts, never sample
names or values.

`postOA.sh` has four steps, all writing to `res/prism/db/`. Rule for
all (user decision): do not infer biology; faithfully preserve what is
there. No record is removed and every value a tool wrote can be
recovered exactly; anything added is a documented extra column, and a
confidence tier may only be a column, never a filter. Errors never pass
silently: every assumption is checked on every row and a failed check
stops the run, naming the column.

1. `scripts/postOA.R` stacks each per-sample table of every group into
   `data/<tool>_<mode>_<table>.tsv.gz` with `group_id` and `sample_id`
   in front. It copies lines as the tool wrote them (no parsing, no
   filtering, no reformatting) and stops on a missing file, a header
   that differs between samples, or a line with a different number of
   fields than its header. Each table has a row in `linx_tables`
   (description, source path), from which `manifest.tsv` is written;
   change the description there when a table or its checked behavior
   changes.
2. `scripts/svCohort.R` writes `sv_cohort_filtered.tsv.gz`, the core
   analysis table: one row per event of PURPLE's somatic
   `purple.sv.vcf.gz` (all FILTER values), with every VCF INFO/FORMAT
   field, the ESVEE `breakend.tsv` columns (joined on VCF ID = ESVEE
   `Id`), the LINX tables of step 1 (joined on `vcfIdStart`;
   per-transcript, fusion, link and driver rows as lists), region
   tracks from `scripts/fetchTracks.sh`, the same breakends in other
   groups' PURPLE somatic, PURPLE germline and ESVEE unfiltered calls,
   BAM read counts (`bamEvidence.sh`, `--bam-dir`) and population SV
   sites (`scripts/popSvSites.sh`). It adds `{start,end}_gene_role`
   (driver gene panel `likelihoodType` of each `linx_bnd_gene`
   element) and `fusion_gene_5p`/`_3p` (breakend gene of each fusion's
   5'/3' breakend ID; LINX names enhancer partners `@IGH` etc., so the
   name is checked as `5p_3p` or `@region_3p`). These carry no
   `linx_bnd_`/`linx_fusion_` prefix: `checkSvCohort.R` rebuilds the
   LINX rows from the columns with those prefixes. The table starts
   with a lead block of key columns (`lead_cols`, user decision; see
   `docs/prism_methods.md` Layout). The column dictionary must equal
   `docs/sv_cohort_filtered.columns.csv` (`--dict-ref`), or the run
   stops before writing the table and writes
   `sv_cohort_filtered.columns.new.csv` to review and copy over it.
   The section "unpack, deduplicate, event-level columns"
   reshapes it for analysis (user decisions), each rule checked on
   every row: exact duplicate columns dropped (ESVEE columns repeating
   VCF fields, StrandBias = FORMAT SB to two decimals, MATEID,
   `vcf_normal_sample`); `PURPLE_AF`/`CN`/`CN_CHANGE` reduced to each
   breakend's own value; fixed pairs split (AD ref/alt, CIPOS and
   IHOMPOS lower/upper, TPR five/three prime span); lists that always
   hold one value kept as scalars (a row with more stops the run);
   fields equal on both breakends stored once as unprefixed event-level
   columns. Lists of 1 to N values stay lists. VCF `.` is stored as
   empty; the TSV is written without quoting. It writes
   `sv_cohort_filtered.columns.tsv` (description, type, `list_sep`,
   `list_count_column`, non-empty count; types from `column_type()`)
   and adds the table, its Parquet copy and the dictionary to
   `manifest.tsv`. GenomicRanges is called through its namespace;
   attaching it masks `fs::path` and dplyr verbs. A new column needs a
   branch in `describe()` and, if its type is not inferred correctly,
   in `column_type()`; the run stops on a column without a description.
3. `scripts/tsvToParquet.py` (python3, pyarrow) writes
   `sv_cohort_filtered.parquet`, typed by the dictionary: int32
   integers, lists as arrays, list lengths checked against the count
   column. Any value that does not fit its type stops it (exit 1).
4. `docs/scripts/checkSvCohort.R` undoes every unpacking rule, rebuilds
   every VCF record and ESVEE row, compares them with the source files,
   recomputes the gene roles and fusion genes,
   and checks every Parquet cell against the TSV (`parquet_equal`, with
   nanoparquet). It exits 1 on any failure, which fails the job.

`postOA.sh` checks before it starts that the R package `nanoparquet` is
installed and that `python3` can import `pyarrow` (installed per user:
`python3 -m pip install --user pyarrow`). `bamEvidence.sh` checks that
`python3` can import `pysam`.

Inputs of step 2 made by other jobs: `bamEvidence.sh` (array job, one
task per BAM; `scripts/bamEvidence.py` counts each BAM at every
breakend of the cohort into `res/prism/db/bam_evidence/`; it needs
`samples.tsv` from step 1) and `scripts/popSvSites.sh` (bcftools query
of the population SV VCFs into `TRACK_DIR/<source>/sites.tsv.gz`).

`docs/prism_methods.md` is the methods record of every post-pipeline
step: exact definitions taken from the code, for the client and for
papers. When a post-pipeline script changes behavior, update its
section there, and rerun `docs/scripts/checkSvCohort.R` after a change
to `svCohort.R`. `checkSvCohort.R` does not recompute the computed
columns (region, `xs_*`, `bam_*`, population); after a change to them,
also run `docs/scripts/spotCheckSvCohort.R` as a job and record the
result in `prism_methods.md` section 1.5.

`scripts/svReport.R` turns the SV table into the workbook people read
(`res/prism/sv/prism_sv_report.xlsx`): PASS events only, gene and
fusion columns from LINX only (no other annotation; user decision),
five evidence flags, high confidence = no flag (user decision
2026-11-23), gene and fusion sheets counted in tumors, a `rule` sheet
and a `glossary` sheet. Same-gene LINX fusions (`gene_5p` = `gene_3p`)
join two parts of one gene (the user treats them as cryptic splice
products; LINX calls no splicing), left out of the fusion sheet
(`fusions_real`) and the
gene fusion counts, kept on the events sheets; `genes`
keeps genes with 2 or more tumors in at least one count (user
decisions). It stops if an SV table column it reads is missing, if a
column has no glossary entry or an entry matches no column (all its
checks: `docs/SV_REPORT.md` section 9; what "no flag" does and does
not prove: section 6.7). `docs/SV_REPORT.md` describes it in full,
including the LINX rules behind each gene column; change the two
together, and rerun the report after every `postOA.sh` build. The SV
table is the evidence store behind the report, not the deliverable.

`scripts/orangeSvReport.R` writes `res/prism/sv/orange_sv_report.xlsx`:
the SV records of every ORANGE JSON (`linx.fusions`,
`linx.somaticBreakends`, `linx.germlineBreakends`), ORANGE fields
unchanged, each placed on its SV table event (sample and LINX `svId`;
fusions through their reported `fusion.tsv` row) with the event's
FILTER and SV report flags added. ORANGE 5.0.1 lists reported LINX
fusions and REPORTED breakends and tests no FILTER value; the script
stops if the JSON is not exactly those LINX rows. Run it after
`svReport.R`; described in `docs/SV_REPORT.md` section 13.

Any statement about what a Hartwig tool does goes into
`docs/ACTUAL_METHODS.md` (or `SV_REPORT.md`) only after it is checked
in the tool's source at the tag the run used (sparse clone of
`hartwigmedical/hmftools` at e.g. `linx-v2.3.1`) and, where possible,
against the run's files; cite the class and method. Where the README
disagrees, record it in `ACTUAL_METHODS.md` section 5. The LINX README
is wrong about `svs.tsv` `geneStart`/`geneEnd`, `breakend.tsv` rows and
the artifact rules.

Post-pipeline layout (user decision): deliverables go under
`res/prism/` in the folder that holds `out/` (the project folder), with
the files people read (PDFs, sample lists) at the top of each folder
and machine tables in `data/`. Logs go to `out/logs/other/`, next to
the pipeline's own; nothing else is written into `out/`.
`deliverLocal.sh` is to be redesigned around `res/prism/`; do not
change it until then.

`WORK_DIR` is the Nextflow work folder of an earlier run; pass it to
resume that run. Resume from the same project folder and with the same
samplesheet.

Checks after a change (there is no test suite):

```
bash -n Prism/runOncoanalyser.sh Prism/deliverLocal.sh Prism/qcCohort.sh Prism/cohortReport.sh Prism/postOA.sh Prism/bamEvidence.sh Prism/scripts/fetchTracks.sh Prism/scripts/popSvSites.sh Prism/docs/scripts/countVcfFilters.sh
shellcheck -x Prism/runOncoanalyser.sh Prism/00.SETUP.sh Prism/deliverLocal.sh Prism/qcCohort.sh Prism/cohortReport.sh Prism/postOA.sh Prism/bamEvidence.sh Prism/scripts/fetchTracks.sh Prism/scripts/popSvSites.sh Prism/docs/scripts/countVcfFilters.sh
Rscript --vanilla -e "for (f in c('Prism/scripts/cohortReport.R', 'Prism/scripts/postOA.R', 'Prism/scripts/svCohort.R', 'Prism/docs/scripts/checkLinx.R', 'Prism/docs/scripts/checkSvCohort.R', 'Prism/docs/scripts/spotCheckSvCohort.R', 'Prism/scripts/svReport.R', 'Prism/scripts/orangeSvReport.R', 'Prism/examples/readSvCohort.R')) invisible(parse(f))"
Rscript --vanilla -e "for (f in Sys.glob('Prism/docs/dbScheme/scripts/*.R')) invisible(parse(f))"
python3 -c 'import ast,sys; ast.parse(open(sys.argv[1]).read())' Prism/scripts/tsvToParquet.py
python3 -c 'import ast,sys; ast.parse(open(sys.argv[1]).read())' Prism/scripts/bamEvidence.py
Prism/bin/nextflow lint Prism/conf/iris.config
Prism/bin/nextflow -c Prism/conf/iris.config config -flat
```

`-c` goes before `config`. In the `config -flat` output,
`beforeScript` and `singularity.runOptions` must show `$SLURM_JOB_ID`
literally. `shellcheck` reports SC1091 for `source "$PDIR/conf/ENV.sh"`
(in `runOncoanalyser.sh` and `qcCohort.sh`); that is expected. `lint`
and `config` write `.nextflow.log` in the current folder, which
rotates a run's log there (`.nextflow.log` becomes `.nextflow.log.1`).
Add `-log <file>` before the command (for example
`Prism/bin/nextflow -log /localscratch/soccin/<task>/lint.log lint ...`)
to keep the run logs in place.

To check a resource request against partition and QOS limits without
submitting anything: `/usr/bin/sbatch --test-only -p <partitions>
--qos=<qos> -c <cpus> --mem=<MB>M -t <time> --wrap=true`.

Do not start a pipeline run or submit a job unless the user asks.

## How the pieces fit

`runOncoanalyser.sh` is the only entry point to the pipeline
(`deliverLocal.sh` and `qcCohort.sh` work on its finished output). It:

1. Finds the Prism folder (`PDIR`). `sbatch` runs a copy of the script
   from the SLURM spool folder, so `$0` is wrong inside a job. The
   script tries, in order, the folder of `$0`, `SBATCH_SCRIPT_DIR` (set
   by the user's `~/bin/sbatch` wrapper), and the `Command=` path from
   `scontrol show job`. A candidate is accepted only if it holds
   `conf/iris.config`. If that file is renamed, update `is_prism_dir`
   (copied in `qcCohort.sh`, `cohortReport.sh` and `postOA.sh`).
2. Sources `conf/ENV.sh`, which sets `NXF_VER`, `NXF_HOME`, a new
   timestamped `NXF_WORK`, the Singularity cache and `NXF_OPTS`.
3. Makes a temp folder with `mktemp` under
   `/localscratch/soccin/oncoanalyser` and points `TMPDIR`,
   `SINGULARITY_TMPDIR`, `APPTAINER_TMPDIR` and `-Djava.io.tmpdir` at it.
4. Runs `$PDIR/bin/nextflow run nf-core/oncoanalyser -revision 3.0.0`
   with `conf/genome_GRCh38_hmf.config`, `conf/iris.config` and
   `conf/params.yaml`, always with `-resume`.

Each config file has one job:

- `conf/ENV.sh`: environment variables only. `00.SETUP.sh` also sources
  it, so `NXF_VER` is defined once.
- `conf/iris.config`: infrastructure only (executor, partition choice,
  resource limits, retries, Singularity, task temp folders). No
  pipeline parameters.
- `conf/genome_GRCh38_hmf.config`: local reference data paths.
- `conf/params.yaml`: `mode`, `genome`, `outdir`. No `input`; the
  samplesheet is always a command-line argument.

Values that must stay in agreement across files:

- The Singularity cache path is in `conf/ENV.sh` (three variables) and
  in `singularity.cacheDir` in `conf/iris.config`.
- The temp root `/localscratch/soccin/oncoanalyser` is in
  `runOncoanalyser.sh`, `00.SETUP.sh`, `qcCohort.sh`, `cohortReport.sh`,
  `postOA.sh`, and twice in `conf/iris.config` (`beforeScript` and
  `runOptions`).
- The pipeline version 3.0.0 is in `-revision` in `runOncoanalyser.sh`
  and in the reference path in `conf/genome_GRCh38_hmf.config`.
- `qcCohort.sh` uses the QSEE image of oncoanalyser 3.0.0
  (`hmftools-qsee:1.0--hdfd78af_0`, from `modules/local/qsee/main.nf`)
  and the QSEE 1.0 input file list. Check both on a pipeline upgrade.
- In `conf/iris.config`, `clusterOptions` (QOS) depends on the `queue`
  rule: it drops `--qos=priority` when the queue includes `cpushort`.
  Change them together.
- When a script, config or path changes, update `README.md` and
  `docs/oncoanalyser_3.0.0_WGS_instructions.md` to match.

`scripts/makeOncoanalyserInput.R` builds the samplesheet from
`mapping.csv` (sample to IGO FASTQ folder) and `manifest.csv` (sample
to subject, group and type). It reads lane, library ID and flowcell
from the IGO path layout
`<run>_<flowcell>/Project_*/Sample_*_IGO_*/*_L00N_R1_001.fastq.gz` and
writes one row per FASTQ pair. The flowcell is in the `info` column
because lane numbers repeat across flowcells and oncoanalyser builds
read group IDs from sample, library, lane and flowcell.

## Rules specific to this repo

- No write to a node's `/tmp`, by any script, job or config. The head
  process uses the `mktemp` folder above. Each task gets
  `/localscratch/soccin/oncoanalyser/$SLURM_JOB_ID` from `beforeScript`,
  mounted over the container's `/tmp` and at its own path. Do not set
  `process.scratch = true`; it creates its folder under `/tmp` before
  `beforeScript` runs.
- No automatic cleanup of temp or work folders. The user removes them.
- The head script must not depend on `SLURM_*` variables, because it is
  also run directly. Use `mktemp`, not `$SLURM_JOB_ID`, for its temp
  folder.
- Nothing under `~`. `NXF_HOME`, `NXF_WORK` and the image cache are
  under `/scratch/core001/bic/socci` (the folder is `socci`; the
  username is `soccin`).
- Call `bin/nextflow`, never a global `nextflow`. It is installed by
  `00.SETUP.sh` and is not tracked in git. The version is pinned
  because Nextflow 26.04 made strict syntax the default.
- Do not set `cpus` or `memory` in `executor {}`; since Nextflow 26.10
  they cap SLURM tasks.
- Keep every `#SBATCH` line in the top comment block of a run script.
  The `sbatch` wrapper inserts its exports before the first line that
  is not blank or a comment. The wrapper's temp copies match `.*.sh.*`
  and are ignored by git.
- Do not add `cd "$SLURM_SUBMIT_DIR"` and do not pass `--export=` with
  a value other than `ALL`.
- `errorStrategy` retries every failure once; after that it retries
  every exit code except 1, because the WiGiTS tools exit 1 on real
  failures. Java `OutOfMemoryError` also exits 1, hence the one retry
  (memory grows with `task.attempt`). Check a new pipeline's exit codes
  before reusing this.
- Tasks on `cmobic_cpu` and `cmobic_short` run with QOS `priority`.
  `cpushort` allows only `normal`, and SLURM has
  `EnforcePartLimits=ALL`: a job listed for several partitions must be
  valid in each one. `resourceLimits` (52 CPUs, 975 GB) follows from the
  same rule; see `docs/NF_SCRIPTING_INFO.md` section 6.
- Before resuming a run, change only process directives (`cpus`,
  `memory`, `time`, `queue`, `clusterOptions`, `errorStrategy`); they
  are not part of the task hash. `ext.*` values that a module script
  uses are, so changing them reruns that process and everything after
  it. `-resume` restarts `task.attempt` at 1, so a task that needed a
  retry for memory needs a larger starting value.
- The pipeline name is spelled "oncoanalyser". All other text uses US
  spelling.
- The repo layout follows https://github.com/soccin/adagio: `00.SETUP.sh`,
  `bin/`, `conf/`, `scripts/`, `docs/`.
