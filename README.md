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

   It hard links each group's ORANGE PDF, QSEE QC PDF, LINX HTML
   report, PURPLE circos plot and CUPPA plot from `OUT_DIR` (default `out`) into
   `RESULTS_DIR/PDF`, `RESULTS_DIR/HTML` and `RESULTS_DIR/PNG` (default
   `results`). Both folders must be on the same file system. Rerun it
   to add new files; existing links are skipped.

Outputs go to the project folder, not to `Prism/`:

- `out/`: pipeline results (`outdir` in `conf/params.yaml`)
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
  00.SETUP.sh           one-time setup: installs Nextflow into bin/
  bin/                  nextflow (installed by 00.SETUP.sh, not tracked)
  conf/
    ENV.sh                   Nextflow environment: version, home, work, caches
    iris.config              SLURM and Singularity settings for IRIS
    genome_GRCh38_hmf.config local reference data paths
    params.yaml              pipeline parameters (mode, genome, outdir)
  scripts/
    makeOncoanalyserInput.R  builds the samplesheet
  examples/             templates for mapping.csv and manifest.csv
  docs/                 pipeline instructions and notes
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

No job writes to a node's `/tmp`. The temp folders on `/localscratch`
are not removed automatically.

## Docs

- `docs/oncoanalyser_3.0.0_WGS_instructions.md`: how the pipeline is
  run in WGS mode and why the settings were chosen
- `docs/NF_SCRIPTING_INFO.md`: rules for writing Nextflow run scripts
  and configs on IRIS, and lessons from building Prism. Copy it into
  the next pipeline's wrapper repo.
