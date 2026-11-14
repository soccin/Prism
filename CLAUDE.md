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
to each rule. It uses the paths of the folder Prism was copied from;
the table at its top maps them to Prism paths.

## Commands

Run from the project folder (the parent of `Prism/`):

```
bash Prism/00.SETUP.sh                                  # once per clone: install Nextflow in Prism/bin
Rscript Prism/scripts/makeOncoanalyserInput.R mapping.csv manifest.csv oncoanalyser_input.csv
sbatch ./Prism/runOncoanalyser.sh INPUT.csv [WORK_DIR]  # as a SLURM job
./Prism/runOncoanalyser.sh INPUT.csv [WORK_DIR]         # in the current shell
```

`WORK_DIR` is the Nextflow work folder of an earlier run; pass it to
resume that run. Resume from the same project folder and with the same
samplesheet.

Checks after a change (there is no test suite):

```
bash -n Prism/runOncoanalyser.sh
shellcheck -x Prism/runOncoanalyser.sh Prism/00.SETUP.sh
Prism/bin/nextflow lint Prism/conf/iris.config
Prism/bin/nextflow -c Prism/conf/iris.config config -flat
```

`-c` goes before `config`. In the `config -flat` output,
`beforeScript` and `singularity.runOptions` must show `$SLURM_JOB_ID`
literally. `shellcheck` reports SC1091 for `source "$PDIR/conf/ENV.sh"`;
that is expected.

Do not start a pipeline run or submit a job unless the user asks.

## How the pieces fit

`runOncoanalyser.sh` is the only entry point. It:

1. Finds the Prism folder (`PDIR`). `sbatch` runs a copy of the script
   from the SLURM spool folder, so `$0` is wrong inside a job. The
   script tries, in order, the folder of `$0`, `SBATCH_SCRIPT_DIR` (set
   by the user's `~/bin/sbatch` wrapper), and the `Command=` path from
   `scontrol show job`. A candidate is accepted only if it holds
   `conf/iris.config`. If that file is renamed, update `is_prism_dir`.
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
  `runOncoanalyser.sh`, `00.SETUP.sh`, and twice in `conf/iris.config`
  (`beforeScript` and `runOptions`).
- The pipeline version 3.0.0 is in `-revision` in `runOncoanalyser.sh`
  and in the reference path in `conf/genome_GRCh38_hmf.config`.
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
- `errorStrategy` retries every exit code except 1, because the WiGiTS
  tools exit 1 on real failures. Check a new pipeline's exit codes
  before reusing this.
- The pipeline name is spelled "oncoanalyser". All other text uses US
  spelling.
- The repo layout follows https://github.com/soccin/adagio: `00.SETUP.sh`,
  `bin/`, `conf/`, `scripts/`, `docs/`.
