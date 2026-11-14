# Notes for agents: Nextflow / nf-core run scripts on IRIS

> Copied into Prism on 2026-11-14 from the project folder where the
> pipeline was first set up. The text below uses that folder's paths.
> In Prism the files are:
>
> | Path in this document | Path in Prism |
> |---|---|
> | `runOncoanalyser.sh` | `runOncoanalyser.sh` |
> | `config/iris.config` | `conf/iris.config` |
> | `config/genome_GRCh38_hmf.config` | `conf/genome_GRCh38_hmf.config` |
> | `ENV.sh` | `conf/ENV.sh` |
> | `params.yaml` | `conf/params.yaml` |
> | `./bin/nextflow` | `bin/nextflow`, installed by `00.SETUP.sh` |
>
> The run script now finds these files from the Prism folder, not from
> the current directory. See `README.md`.

Written 2026-11-13, after setting up nf-core/oncoanalyser 3.0.0 in this
folder. Every item below either broke or had to be corrected by the user
in that session.

Before writing a run script or Nextflow config, read this file and both
CLAUDE.md files:
- `~/.claude/CLAUDE.md`
- `/data1/core001/work/bic/socci/CLAUDE.md` (the same file as
  `~/Work/CLAUDE.md`)

Working examples in this folder:
- `runOncoanalyser.sh`
- `config/iris.config`
- `ENV.sh` (the user's file)
- `params.yaml`
- `config/genome_GRCh38_hmf.config`

## 1. Never write to a node's /tmp

This is the most important rule. It covers every script, job and config
you write, not only the commands you run yourself. A write to a node's
`/tmp` damages the node for every user of it.

### What went wrong
The first test run wrote temporary files to the nodes' `/tmp`. The
oncoanalyser code sets no temp location; a grep for `scratch`, `tmpdir`,
`runOptions` and `beforeScript` found nothing. All of the `/tmp` use came
from defaults:
- **Head process:**
  - The Nextflow JVM: `java.io.tmpdir` defaults to `/tmp`. This holds the
    file collector folders `nxf-*` and plugin downloads.
  - Singularity image pulls.
- **Tasks:** `/etc/apptainer/apptainer.conf` has `mount tmp = yes`, so
  the node's `/tmp` is mounted into every container. The tools inside
  (the Java tools, sambamba, STAR, ...) write there by default.

Assume every new pipeline needs the setup below.

### Head process (in the run script)
```bash
mkdir -p /localscratch/soccin/<pipeline>
NXF_LOCAL_TMP=$(mktemp -d -p /localscratch/soccin/<pipeline> head_XXXXXX)
export TMPDIR=$NXF_LOCAL_TMP
export SINGULARITY_TMPDIR=$NXF_LOCAL_TMP
export APPTAINER_TMPDIR=$NXF_LOCAL_TMP
export NXF_OPTS="${NXF_OPTS:-} -Djava.io.tmpdir=$NXF_LOCAL_TMP"
```
- **Use `mktemp`, not `$SLURM_JOB_ID`.** The user runs the script with
  `sbatch` and also directly (`./runX.sh`) for tests.
  - Run directly inside the interactive allocation, `SLURM_JOB_ID` is
    that session's ID: 8901522, the PLEIADES session on isca006. Every
    run would get that same ID.
  - Outside any allocation, `SLURM_JOB_ID` is unset, and `set -u` stops
    the script.
  - The first version used `$SLURM_JOB_ID` and had to be redone.
- **Java ignores `TMPDIR`** and needs `-Djava.io.tmpdir`.
  - `./bin/nextflow` sets `NXF_TEMP=${NXF_TEMP:-$TMPDIR}` and passes
    `-Djava.io.tmpdir=$NXF_TEMP`, so `TMPDIR` already covers the launcher.
  - Set `-Djava.io.tmpdir` explicitly anyway.

### Task jobs (in the Nextflow config)
```groovy
process {
    beforeScript = '''
        NXF_LOCAL_TMP=/localscratch/soccin/<pipeline>/$SLURM_JOB_ID
        mkdir -p $NXF_LOCAL_TMP
        export TMPDIR=$NXF_LOCAL_TMP
        export SINGULARITY_TMPDIR=$NXF_LOCAL_TMP
        export APPTAINER_TMPDIR=$NXF_LOCAL_TMP
    '''
}
singularity {
    runOptions = '-B /localscratch/soccin/<pipeline>/$SLURM_JOB_ID:/tmp -B /localscratch/soccin/<pipeline>/$SLURM_JOB_ID'
}
```
- **Job ID:** each task is its own SLURM job, so `$SLURM_JOB_ID` is
  always set in a task and is unique to it.
- **Single quotes:** Groovy leaves `$SLURM_JOB_ID` alone, and the shell on
  the node fills it in. Check this with `config -flat` (section 4).
- **Order:** `beforeScript` runs on the node in `nxf_main`, before the
  container starts, so the folder exists when it is mounted.
- **The exports are required.**
  - Nextflow submits tasks with sbatch's default `--export=ALL`, so each
    task inherits the head's `TMPDIR`. That folder is on the head's node
    and does not exist on the task's node.
  - `.command.run` passes `TMPDIR` (and `TMP`, if set) into the container
    as `SINGULARITYENV_TMPDIR`. The exports replace the inherited value.
  - The first version had only the `mkdir` and had to be redone.
- **Two mounts:**
  - The first, `-B ...:/tmp`, catches tools that write to `/tmp`
    directly and ignore `TMPDIR`.
  - The second mounts the folder at its own path, so the `TMPDIR` path
    exists inside the container. `/localscratch` is not mounted into
    containers otherwise.
- **`TMP`** was not set in this environment (checked 2026-11-13). If it
  ever is, export it as well.
- **Do not set `process.scratch = true`.**
  - In `.command.run`, the scratch folder is made with
    `nxf_mktemp $TMPDIR`, which defaults to `/tmp`. This happens before
    `beforeScript` runs.
  - If node-local scratch is wanted, give an explicit path under
    `/localscratch/soccin/<pipeline>`.
  - This is read from the `.command.run` template and has not been
    tested.

### Do not
- Never set `TMPDIR=/tmp`, inside or outside containers.
- No automatic cleanup, and do not offer it:
  - no `afterScript 'rm -rf ...'`
  - no `rm -rf` at the end of a script

  The user removes the folders by hand.
- `/localscratch` is node-local, and nodes cannot see each other's.
  - Use it for temporary files only.
  - Anything another job reads goes on shared storage
    (`/scratch/core001/bic/socci/`, `/data1`).

### Check on the first test run
1. `./bin/nextflow -c config/iris.config config -flat`: `beforeScript`
   and `runOptions` show `$SLURM_JOB_ID` literally.
2. The run prints `TMPDIR=/localscratch/soccin/<pipeline>/head_...`.
3. In `.nextflow.log`, `grep localscratch` shows the `nxf-*` collector
   folders there.
4. In the first task's `.command.run`:
   - The `beforeScript` block (`mkdir` and the 3 exports) runs before
     `nxf_launch`.
   - The `singularity exec` line has both `-B` mounts.

Status on 2026-11-13: checks 1 to 4 passed, on work folder
`/scratch/core001/bic/socci/Oncoanalyser/20261113_131345`. Not checked:
what is actually in a task node's `/localscratch`. That needs a shell on
that node, and `srun` is not allowed.

### The TMPDIR rule mix-up
- An older CLAUDE.md line said "Do not set TMPDIR or CLAUDE_CODE_TMPDIR
  to /tmp". It was about the agent's own shell.
- The agent made two mistakes with it:
  - cited it as a reason not to set `TMPDIR` in the pipeline jobs
  - then proposed `TMPDIR=/tmp` inside containers
- Both were wrong. Jobs must set `TMPDIR` to `/localscratch`.
- The current CLAUDE.md says the rule also covers scripts and configs.

## 2. No home folder; use the user's paths

- **Nothing under `~`.** The user had to repeat this rule. Defaults that
  write there:
  - `NXF_HOME` (`~/.nextflow`), which holds the Nextflow runtime, the
    pipeline code and plugins.
  - `SINGULARITY_CACHEDIR` and `APPTAINER_CACHEDIR`.
- **The project's `ENV.sh`** (the user's file) sets:
  ```bash
  export NXF_VER=25.10.4
  export NXF_HOME=/scratch/core001/bic/socci/Oncoanalyser/nextflow_home
  export NXF_WORK=/scratch/core001/bic/socci/Oncoanalyser/$(date +%Y%m%d_%H%M%S)
  export NXF_SINGULARITY_CACHEDIR=/scratch/core001/bic/socci/singularity-cache
  export SINGULARITY_CACHEDIR=/scratch/core001/bic/socci/singularity-cache
  export APPTAINER_CACHEDIR=/scratch/core001/bic/socci/singularity-cache
  export NXF_OPTS='-Xms1g -Xmx4g'
  ```
- **nf-core work** goes in `/scratch/core001/bic/socci/<Pipeline>/`.
  - Capitalize the folder name to match the existing ones (`SMap`,
    `Piano`, `Adagio`, `Oncoanalyser`).
  - `NXF_HOME` goes in the same folder.
- **Container images:** reuse `/scratch/core001/bic/socci/singularity-cache`.
  - Do not create a second cache.
  - `singularity.cacheDir` in the config must match.
- **oncoanalyser reference data:**
  `/scratch/core002/nds/nf-core/reference_data/3.0.0/`, downloaded with
  `--outdir /scratch/core002/nds/nf-core`.
- **Do not invent paths.** Wrong paths suggested in the session:
  - `/scratch/core001/bic/soccin/...`: the folder is `socci`, even though
    the username is `soccin`.
  - `/scratch/core002/nds/nf-core/work`.

  Run `ls` to see what exists, or ask.
- **Launch from the project folder.** `.nextflow/` and `.nextflow.log`
  are written to the folder you launch from.

## 3. Nextflow version

- **Always use `./bin/nextflow`,** never a global `nextflow`.
  - The user installs a pinned launcher for each project or pipeline: in
    `bin/`, run `export NXF_VER=25.10.4`, then
    `curl -s https://get.nextflow.io | bash`.
  - The launcher then defaults to that version
    (`NXF_VER=${NXF_VER:-'25.10.4'}`).
- **Why it is pinned:** Nextflow 26.04 made strict syntax the default,
  which broke older scripts.
  - Check the pipeline's `manifest.nextflowVersion`.
  - oncoanalyser 3.0.0 has `'!>=25.10.4'`, with no maximum.
- **Nextflow 26.10:** `executor.cpus` and `executor.memory` also cap SLURM
  jobs. Do not set them in `executor {}`.

## 4. Checking configs and scripts

- **Syntax:** `./bin/nextflow lint <file>`. The user said "lint"; do not
  read it as "link".
- **Resolved values:** `./bin/nextflow -c <file> config -flat`. `-c` goes
  before `config`; `nextflow config -c <file>` is wrong.
- **If `NXF_HOME` does not exist** (the user may have deleted it), do not
  recreate the user's folder. Use a throwaway one:
  `NXF_VER=25.10.4 NXF_HOME=/localscratch/soccin/<task>/nxf_home ./bin/nextflow lint <file>`.
  Both commands download the runtime into `NXF_HOME` the first time.
- **Shell scripts:** `bash -n <script>`.
- **"Write X, do not run"** means: write it, check it as above, and stop.
  Do not execute it or submit it.

## 5. SLURM and sbatch

### Submitting
- Submit with plain `sbatch script.sh ARGS`.
- Do not put `--parsable` in run instructions. It only changes the output
  to the bare job ID. Use it only when a script captures the ID, as in
  the CLAUDE.md pattern for waiting on a job.
- Do not add `cd "$SLURM_SUBMIT_DIR"`.
  - Batch jobs already start in the submit folder.
  - If the user submits with `-D <dir>`, the `cd` undoes that choice.

### The user's wrapper, `~/bin/sbatch`
The `sbatch` on PATH is this wrapper. It:
- copies the script to a hidden temp file next to it
- inserts `export SBATCH_SCRIPT_DIR=<script dir>` and
  `export IRIS_VAR1=...` before the first line that is not blank or a
  comment
- calls `/usr/bin/sbatch`

What this means for scripts:
- Keep all `#SBATCH` lines in the top comment block.
- Use `$SBATCH_SCRIPT_DIR` to find files next to the script. `$0` inside
  the job is SLURM's copy of the script.
- Do not pass `--export=` with anything other than `ALL`. Since April 2026
  that triggers a broken user-environment probe on isca001 / bic_devs.
- Nextflow's task submissions also go through this wrapper. This works:
  `scontrol write batch_script` for task 10652310 showed the inserted
  exports.

### Where jobs run
- **Never use `srun`.** Sessions run inside an existing allocation, so
  `srun` starts a step in that allocation (often on isca006) instead of a
  new job.
- **isca006** (bic_devs, interactive): run nothing heavy there; use
  `sbatch`.
  - The user may choose to run a short test head process there directly
    while watching the load. That is their decision.
  - State the rule once and accept their answer.

### Partitions and limits
- **Walltime limits:** cpushort 2h, cmobic_short 3h, cmobic_cpu 7 days.
- **Choose by the time each task requests,** as in `config/iris.config`:
  ```groovy
  queue = {
      if( !task.time )
          return 'cmobic_cpu'
      if( task.time <= 2.h )
          return 'cpushort,cmobic_short'
      if( task.time <= 3.h )
          return 'cmobic_short'
      return 'cmobic_cpu'
  }
  resourceLimits = [ cpus: 32, memory: 120.GB, time: 7.d ]
  ```
- **The time cap is 7 days.** An earlier draft had 48h.
- **Head job header:** `-p cmobic_cpu -c 2 --mem=8G -t 7-00:00:00`.

## 6. Run script conventions (see runOncoanalyser.sh)

- **Arguments:**
  - Explicit, required arguments: the input samplesheet is the first
    argument. No argument-free commands.
  - Keep `input` out of `params.yaml`. A command-line `--input` overrides
    the params file anyway.
- **Help and errors:**
  - `-h` / `--help` prints usage and exits before anything else runs.
  - A wrong argument count prints usage to stderr and exits 1.
- **Checks:** `set -eu`. Check that the input file exists (`realpath`),
  and that the optional `WORK_DIR` exists.
- **It must work both via `sbatch` and run directly.** Nothing in the
  head script may depend on `SLURM_*` variables.
- **Paths:** `source ENV.sh`; do not copy its exports into the script.
- **Output:** print `INPUT`, `NXF_WORK` and `TMPDIR` at the start.
- **Finish:** `chmod +x`, then check with `bash -n`.
- **Usage line:** `sbatch runOncoanalyser.sh INPUT.csv [WORK_DIR]`.

## 7. Work folders and -resume

- **A new work folder per run:** `NXF_WORK` gets a new timestamped folder
  each time `ENV.sh` is sourced.
- **`-resume` only reuses tasks from the same folder,** so the script
  takes the old folder as an optional second argument. Also:
  - Resume with the same samplesheet.
  - `-resume` needs `.nextflow/` (the run history) in the launch folder.
    Without it, Nextflow prints "Option `-resume` is ignored".
- **Large remote files** go into the work folder and can be cut off.
  - Example: the 27 GB STAR index tarball. Nextflow retries the whole
    file.
  - Download big references yourself instead, with `curl -C -` in an
    sbatch job.

## 8. Config layout

- `config/iris.config`: infrastructure only (executor, queue, limits,
  retries, Singularity, temp setup). No pipeline parameters.
- `config/genome_<genome>.config`: reference paths.
- `params.yaml`: `mode`, `genome`, `outdir`. No `input`.
- `ENV.sh`: environment variables (the user's file).
- **oncoanalyser `errorStrategy`:** its tools exit 1 on real failures,
  hence `{ task.exitStatus != 1 ? 'retry' : 'finish' }` with
  `maxRetries = 3`. Check which exit codes a new pipeline's tools use
  before copying this.

## 9. Working with this user

- **Files:**
  - Do not edit files you did not create (CLAUDE.md, `ENV.sh`,
    `~/bin/sbatch`) unless asked. Propose the text instead.
  - Never delete anything you did not create, and never use `rm -r`.
- **Commands:** one command per Bash call; no `&&`, `;` or `|`.
- **Replies:**
  - Keep them short.
  - When the user is confused, give a numbered plan of a few lines, not
    more explanation.
  - Proposed CLAUDE.md text must be short, because it costs context in
    every session.
- **Scope:**
  - Do not add things nobody asked for. Automatic cleanup was rejected
    twice.
  - Once the user decides something, apply it; do not reopen it.
- **Report:** after every script or config change, update
  `oncoanalyser_3.0.0_WGS_instructions.md` to match.
- **Spelling:** the pipeline is "oncoanalyser" (a proper name); use US
  spelling for everything else. No emoji.
- **Certainty:** ask until you are 95% sure. Do not guess paths,
  partitions or limits.

## 10. Checklist for a new pipeline

1. Read both CLAUDE.md files and this file.
2. Use a pinned `./bin/nextflow`, and check `manifest.nextflowVersion`.
3. In `ENV.sh`, put `NXF_HOME`, `NXF_WORK` and the caches under
   `/scratch/core001/bic/socci/<Pipeline>`. Nothing goes in `~`.
4. Add the head temp setup (section 1).
5. Add the task temp setup: `beforeScript` and the two `-B` mounts
   (section 1).
6. Add the queue closure and `resourceLimits` with time `7.d`.
7. Write the run script as in section 6.
8. Check with `bash -n`, `nextflow lint` and `nextflow -c ... config -flat`.
9. Do not run anything. Tell the user to check the first task's
   `.command.run` (section 1, "Check on the first test run").
