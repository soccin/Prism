# Prism

Wrapper script to run nf-core/oncoanalyser

Runs [nf-core/oncoanalyser](https://nf-co.re/oncoanalyser/3.0.0) 3.0.0
(WGS/WTS mode, genome GRCh38_hmf) on IRIS (SLURM, Singularity).

## Install

Clone Prism into the project folder, then install Nextflow into
`Prism/bin`:

```
cd <project folder>
git clone https://github.com/soccin/Prism.git Prism
./Prism/00.SETUP.sh
```

`00.SETUP.sh` installs the Nextflow version pinned by `NXF_VER` in
`conf/ENV.sh`. `bin/nextflow` is not tracked in git.

## Run

All commands are run from the project folder.

1. Write `mapping.csv` and `manifest.csv` (templates in `examples/`),
   then build the samplesheet:

   ```
   Rscript Prism/scripts/makeOncoanalyserInput.R mapping.csv manifest.csv oncoanalyser_input.csv
   ```

2. Start the pipeline, as a SLURM job or in the current shell:

   ```
   sbatch ./Prism/runOncoanalyser.sh oncoanalyser_input.csv
   ./Prism/runOncoanalyser.sh oncoanalyser_input.csv
   ```

   `./Prism/runOncoanalyser.sh -h` prints the usage.

3. After the run, collect the reports for delivery:

   ```
   ./Prism/deliverLocal.sh [OUT_DIR] [RESULTS_DIR]
   ```

   It hard links each group's ORANGE PDF, QSEE QC report (PDF and
   PNG), LINX HTML report, PURPLE circos plot and CUPPA plot from
   `OUT_DIR` (default `out`) into `RESULTS_DIR/PDF`,
   `RESULTS_DIR/HTML` and `RESULTS_DIR/PNG` (default `results`). Both
   folders must be on the same file system. Rerun it to add new
   files; existing links are skipped.

4. Run QSEE (Hartwig's QC tool) once over all groups, as a SLURM job
   or on a compute node:

   ```
   sbatch ./Prism/qcCohort.sh [OUT_DIR]
   ```

   oncoanalyser runs QSEE once per group. This gives one QC status
   table, one metric table and one PDF (one page per tumor) for the
   whole run, in `res/prism/qc/hmf/` next to `OUT_DIR`: the PDFs and
   sample lists at the top, the tables in `data/`. Its logs go to
   `OUT_DIR/logs/other/`. Groups with and without a normal are run
   separately (`tumor_normal`, `tumor_only`), because QSEE needs a
   normal for every tumor or for none.

5. Write the cohort report, as a SLURM job or on a compute node:

   ```
   sbatch ./Prism/cohortReport.sh [-m MANIFEST] [-g GENES] [OUT_DIR]
   ```

   It reads the result tables of every group and writes
   `res/prism/cohort/prism_report.xlsx`, with one CSV per table in
   `data/`: the sample-level results with their counts over the
   cohort, and the events per sample (somatic variants, gene copy
   number, gene disruptions with breakpoints and support, fusions,
   viruses, whole genome duplication, DPYD status). It counts and
   lists; it applies no significance test and leaves no sample out.
   The first sheet is for the genes of interest; the last two
   are `manifest`, which says what each sheet holds, and `glossary`,
   which explains every column, once.

   - `MANIFEST`: CSV with `sample_id` (tumor sample) and
     `sample_class` (for example `Tumor`, `CellLine`); other columns
     are copied to the `samples` sheet. Without it every sample has
     class `Tumor`.
   - `GENES`: text file, one gene symbol per line: the genes of
     interest. They get the first sheet, `genes_of_interest_events`:
     every event in them, one row per event.
   - Without `-g` the genes of interest are the genes with a pipeline
     driver call of likelihood 0.8 or higher in at least one sample.
     That list is written to `genes_of_interest_default.txt` next to
     the workbook; edit a copy and give it back with `-g`.
   - Events on chrX and chrY are reported only for samples with a
     value in a manifest column named `sex` and a matched normal.
     Events in a gene of `GENES` are reported for every sample.
   - A breakpoint found at the same base pair in more than one sample
     is flagged (`same_bp_flag`) and listed in a workbook of its own,
     `prism_sv_same_breakpoint.xlsx`.

   It needs `Rscript` on the `PATH`, with the R packages `tidyverse`,
   `fs` and `openxlsx`. A rerun replaces the workbooks, the CSVs and,
   without `-g`, the default gene list.

6. Build the cohort tables, as a SLURM job or on a compute node. Once
   per site, download the region tracks first (login node, downloads
   only):

   ```
   ./Prism/scripts/fetchTracks.sh [TRACK_DIR]
   sbatch ./Prism/postOA.sh [-t TRACK_DIR] [OUT_DIR]
   ```

   `fetchTracks.sh` writes GIAB v3.6 GRCh38 stratifications (mappability,
   segmental duplications, tandem repeats, satellites, GC extremes, gaps
   and other difficult regions), the ENCODE blacklist v2 and UCSC
   RepeatMasker and centromeres to `TRACK_DIR`, with
   `TRACK_DIR/tracks.tsv` (source URL, date, md5, GIAB md5 check).

   `postOA.sh` writes to `res/prism/db/`:
   - `data/linx_*.tsv.gz`: each per-sample LINX table of every group
     stacked into one table, with `group_id` and `sample_id` in front.
     Lines are copied as LINX wrote them. Somatic `svs`, `clusters`,
     `links`, `breakend`, `fusion`, `drivers` and `driver_catalog` for
     every group; germline `svs`, `breakend`, `disruption`, `clusters`,
     `links` and `driver_catalog` for groups with a matched normal.
   - `sv_cohort_filtered.tsv.gz`: one row per somatic SV event: every
     record of PURPLE's `purple.sv.vcf.gz` (all FILTER values; nothing
     left out). Each event has every field PURPLE (VCF INFO and FORMAT
     of tumor and normal), ESVEE (`breakend.tsv`) and LINX (SV, cluster,
     breakends per transcript, fusions, links, cluster drivers) wrote
     for it, the region tracks at each breakend, and whether the same
     breakend or pair is in other groups' calls. Values are copied, not
     interpreted.
   - `sv_cohort_filtered.columns.tsv`: the source and description of
     every column of that table, and its number of non-empty values.
   - `manifest.tsv` (what each file holds, where it comes from, its
     groups and rows), `samples.tsv` (tumor, normal, run mode and LINX
     version per group) and `rows_per_sample.tsv`.

   It stops if a file is missing, if the files of a table do not match,
   or if the SV events do not account for every PURPLE record and LINX
   row. It needs `Rscript` with `tidyverse`, `fs` and the Bioconductor
   packages `GenomicRanges`, `IRanges` and `S4Vectors`. A rerun
   replaces the files.

Outputs go to the project folder, not to `Prism/`:

- `out/`: pipeline results (`outdir` in `conf/params.yaml`)
- `res/prism/`: deliverables made after the pipeline from `out/`
  (`res/prism/qc/hmf/` from `qcCohort.sh`, `res/prism/cohort/` from
  `cohortReport.sh`, `res/prism/db/` from `postOA.sh`); their logs
  go to `out/logs/other/`
- `oncoanalyser.<jobid>.out`: head job log (sbatch only)
- `.nextflow.log`, `.nextflow/`: Nextflow log and run history

The Nextflow work folder is a new
`/scratch/core001/bic/socci/Oncoanalyser/<YYYYMMDD_HHMMSS>` for each
run. The log starts with four lines:

- `PRISM=`: the Prism folder and its git version (`git describe`)
- `INPUT=`: the samplesheet, as an absolute path
- `NXF_WORK=`: the work folder
- `TMPDIR=`: the head job temp folder

## Resume

Pass the work folder of the earlier run as the second argument:

```
sbatch ./Prism/runOncoanalyser.sh oncoanalyser_input.csv /scratch/core001/bic/socci/Oncoanalyser/<YYYYMMDD_HHMMSS>
```

Resume from the same project folder as the first run, because Nextflow
keeps its run history in `.nextflow/` there. `-resume` continues the
last run in that history, so do not start another `nextflow run` from
the project folder in between.

If a task failed for lack of memory, raise its starting resources in
`conf/iris.config` before resuming: on resume, retries start again at
attempt 1. Change only resources (`cpus`, `memory`, `time`), the
partition or QOS settings, or `errorStrategy`. These do not change the
task hash, so finished tasks stay cached. Changing an `ext.*` setting
does, and reruns the process and everything after it.

In the first minutes of a resumed run, `.nextflow.log` should list the
finished tasks as `Cached process`; only unfinished tasks should be
`Submitted`.

## Layout

```
Prism/
  runOncoanalyser.sh    run script (sbatch or direct)
  deliverLocal.sh       hard links the reports from out/ into results/
  qcCohort.sh           QSEE over all groups, into res/prism/qc/hmf/
  cohortReport.sh       cohort report, into res/prism/cohort/
  postOA.sh             cohort tables (stacked LINX tables, SV event table), into res/prism/db/
  00.SETUP.sh           one-time setup: installs Nextflow into bin/
  bin/                  nextflow (installed by 00.SETUP.sh, not tracked)
  conf/
    ENV.sh                   Nextflow environment: version, home, work, caches
    iris.config              SLURM and Singularity settings for IRIS
    genome_GRCh38_hmf.config local reference data paths
    params.yaml              pipeline parameters (mode, genome, outdir)
  scripts/
    makeOncoanalyserInput.R  builds the samplesheet
    cohortReport.R           builds the cohort report (run by cohortReport.sh)
    postOA.R                 stacks the LINX tables (run by postOA.sh)
    svCohort.R               builds the SV event table (run by postOA.sh)
    fetchTracks.sh           downloads the region tracks postOA.sh uses
  examples/             templates for mapping.csv and manifest.csv
  docs/                 pipeline instructions and notes
    dbScheme/           which file holds each per-sample value
```

## How the run script finds the Prism folder

`sbatch` runs a copy of the script from the SLURM spool folder, so the
script cannot use its own path there. It tries these in order and uses
the first folder that holds `conf/iris.config`:

1. the folder of the script itself (direct run)
2. `SBATCH_SCRIPT_DIR` (set by the `~/bin/sbatch` wrapper)
3. the path of the submitted script, from `scontrol show job` (plain
   `sbatch`)

## Site settings

The paths are specific to IRIS and to user `soccin`:

- `conf/ENV.sh`: `NXF_HOME`, `NXF_WORK` and the Singularity cache under
  `/scratch/core001/bic/socci`
- `conf/iris.config`: partitions, QOS `priority` on the cmobic
  partitions, resource caps (52 CPUs, 975 GB), more memory or CPUs for
  PURPLE, ESVEE and SAGE somatic, retries, the Singularity cache, and
  temp folders under `/localscratch/soccin/oncoanalyser`
- `conf/genome_GRCh38_hmf.config`: reference data under
  `/scratch/core002/nds/nf-core/reference_data/3.0.0`
- `runOncoanalyser.sh`: partition `cmobic_cpu`, QOS `priority`, head job temp folder
  under `/localscratch/soccin/oncoanalyser`
- `postOA.sh` and `scripts/fetchTracks.sh`: region tracks under
  `/scratch/core001/bic/socci/reference/prism_tracks/GRCh38`

No job writes to a node's `/tmp`. The temp folders on `/localscratch`
are not removed automatically.

## Docs

- `docs/oncoanalyser_3.0.0_WGS_instructions.md`: how the pipeline is
  run in WGS mode and why the settings were chosen
- `docs/NF_SCRIPTING_INFO.md`: rules for writing Nextflow run scripts
  and configs on IRIS, and lessons from building Prism. Copy it into
  the next pipeline's wrapper repo.
- `docs/ACTUAL_METHODS.md`: what the pipeline actually does at each
  step, checked against the code and the run's own logs and files,
  where it differs from the tool READMEs, and how to check it again
  (`docs/scripts/countVcfFilters.sh`, `docs/scripts/checkLinx.R`)
- `docs/LINX_tables.md`: why LINX writes seven somatic SV tables,
  what one row of each is, how they link, and R examples on the
  `res/prism/db` tables
- `docs/OA_AlgoKey.md`: what each Hartwig tool name (PURPLE, AMBER,
  ...) means in plain words
- `docs/dbScheme/`: which output file first holds each per-sample
  value (primary) and which files only copy it (ORANGE, QSEE), field by
  field, as the basis for the cohort tables (`postOA.sh`).
  Covers the per-sample summary files and the ORANGE JSON so far; see
  its `README.md` for the files and how to rebuild them.
