# Notes for agents: wrapping an nf-core pipeline on IRIS

Lessons from building Prism (https://github.com/soccin/Prism), the
wrapper for nf-core/oncoanalyser 3.0.0. First written 2026-11-13, when
the pipeline was set up directly in a project folder. Rewritten
2026-11-14, after that setup was moved into the Prism repo.

Each new pipeline gets its own wrapper repo with its own name, and this
file is copied into that repo's `docs/`. Read all of it, and these,
before writing a run script or a Nextflow config:
- `~/.claude/CLAUDE.md`
- `/data1/core001/work/bic/socci/CLAUDE.md` (the same file as
  `~/Work/CLAUDE.md`)
- https://github.com/soccin/adagio, the user's pattern for wrapper repos
- Prism itself, the working example for every file named below

Where the rules come from: items under "What went wrong" broke or were
corrected by the user. Items marked "Decision" are the user's choices.
The rest describe how Prism does it. Section 14 lists what has not been
tested. Section 16 is the checklist.

Paths and names below are Prism's. Section 11 lists every place that
changes for a new pipeline.

## 1. Repo layout

Decision (2026-11-14):
- One repo per pipeline. The user picks the name.
- The repo is cloned into a project folder. Every command runs from the
  project folder, not from inside the repo. Outputs (`out/`,
  `.nextflow/`, `.nextflow.log`, the job log) go to the project folder.
  Nothing a run writes belongs in the repo.
- The run script must work both under `sbatch` and run directly.
- The layout follows adagio:

```
Prism/
  runOncoanalyser.sh    run script (sbatch or direct)
  00.SETUP.sh           installs Nextflow into bin/, once per clone
  bin/                  nextflow (installed, not tracked)
  conf/
    ENV.sh                   environment variables
    iris.config              SLURM and Singularity settings
    genome_GRCh38_hmf.config reference data paths
    params.yaml              pipeline parameters
  scripts/              helpers (the samplesheet builder)
  examples/             input templates
  docs/                 pipeline instructions and this file
  README.md
  CLAUDE.md
```

Notes:
- Prism keeps the run script at the repo top because the user asked for
  `./Prism/runOncoanalyser.sh INPUT.csv`. adagio keeps its run scripts
  in `bin/`. Ask the user which the new repo uses.
- `.gitignore` has `bin/nextflow` and `.*.sh.*`. The second pattern
  covers the temp copies the `sbatch` wrapper writes next to the script
  (section 6).
- `CLAUDE.md` (made with `/init`) lists the commands, how the files fit
  together, the values that must agree across files, and the repo's
  rules.
- Prism is public. Its `examples/` hold placeholder values only
  (`SUBJECT01`, `Project_00000_A`), and no file names the project the
  pipeline was first run for.

### Git and GitHub
- The session's working folder is usually the project folder, which is
  not a git repo. Use `git -C <repo> ...`.
- The user made the Prism folder and its first commit ("primum") before
  the work started. Run `git log` first.
- Work on a branch (`feat/setup-01`). Commit in logical groups, only
  when asked, with the user's `/draft-commit` skill. Merge to master
  (for Prism this was a fast-forward).
- A push publishes the code. Ask "Are you sure?" every time and wait.
  Prism was created with
  `gh repo create Prism --public --source <repo> --remote origin`, then
  `git push -u origin master`. Ask whether the new repo is public or
  private. Do not guess its name or owner.
- The clone URL for the README exists only after the GitHub repo does.
  Fill it in then.

## 2. Never write to a node's /tmp

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

### Head process (in the run script and in 00.SETUP.sh)
```bash
mkdir -p /localscratch/soccin/<pipeline>
NXF_LOCAL_TMP=$(mktemp -d -p /localscratch/soccin/<pipeline> head_XXXXXX)
export TMPDIR=$NXF_LOCAL_TMP
export SINGULARITY_TMPDIR=$NXF_LOCAL_TMP
export APPTAINER_TMPDIR=$NXF_LOCAL_TMP
export NXF_OPTS="${NXF_OPTS:-} -Djava.io.tmpdir=$NXF_LOCAL_TMP"
```
- **Use `mktemp`, not `$SLURM_JOB_ID`.** The user runs the script with
  `sbatch` and also directly (`./Prism/runX.sh`) for tests.
  - Run directly inside the interactive allocation, `SLURM_JOB_ID` is
    that session's ID: 8901522, the PLEIADES session on isca006. Every
    run would get that same ID.
  - Outside any allocation, `SLURM_JOB_ID` is unset, and `set -u` stops
    the script.
  - The first version used `$SLURM_JOB_ID` and had to be redone.
- **Java ignores `TMPDIR`** and needs `-Djava.io.tmpdir`.
  - `bin/nextflow` sets `NXF_TEMP=${NXF_TEMP:-$TMPDIR}` and passes
    `-Djava.io.tmpdir=$NXF_TEMP`, so `TMPDIR` already covers the launcher.
  - Set `-Djava.io.tmpdir` explicitly anyway.

### Task jobs (in conf/iris.config)
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
  the node fills it in. Check this with `config -flat` (section 5).
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
1. `Prism/bin/nextflow -c Prism/conf/iris.config config -flat`:
   `beforeScript` and `runOptions` show `$SLURM_JOB_ID` literally.
2. The run prints `TMPDIR=/localscratch/soccin/<pipeline>/head_...`.
3. In `.nextflow.log`, `grep localscratch` shows the `nxf-*` collector
   folders there.
4. In the first task's `.command.run`:
   - The `beforeScript` block (`mkdir` and the 3 exports) runs before
     `nxf_launch`.
   - The `singularity exec` line has both `-B` mounts.

For results so far, see section 14.

### The TMPDIR rule mix-up
- An older CLAUDE.md line said "Do not set TMPDIR or CLAUDE_CODE_TMPDIR
  to /tmp". It was about the agent's own shell.
- The agent made two mistakes with it:
  - cited it as a reason not to set `TMPDIR` in the pipeline jobs
  - then proposed `TMPDIR=/tmp` inside containers
- Both were wrong. Jobs must set `TMPDIR` to `/localscratch`.
- The current CLAUDE.md says the rule also covers scripts and configs.

## 3. No home folder; use the user's paths

- **Nothing under `~`.** The user had to repeat this rule. Defaults that
  write there:
  - `NXF_HOME` (`~/.nextflow`), which holds the Nextflow runtime, the
    pipeline code and plugins.
  - `SINGULARITY_CACHEDIR` and `APPTAINER_CACHEDIR`.
- **`conf/ENV.sh`** is a copy of the user's project `ENV.sh`. It sets:
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
- **Do not invent paths.** Wrong paths suggested in the first session:
  - `/scratch/core001/bic/soccin/...`: the folder is `socci`, even though
    the username is `soccin`.
  - `/scratch/core002/nds/nf-core/work`.

  Run `ls` to see what exists, or ask.
- **Launch from the project folder.** `.nextflow/` and `.nextflow.log`
  are written to the folder you launch from.

## 4. Nextflow version

- **`conf/ENV.sh` pins it** with `NXF_VER`; no other file sets the
  version. `00.SETUP.sh` installs that version into `bin/` (section 8).
- **Always call the repo's `bin/nextflow`,** never a global `nextflow`:
  `"$PDIR/bin/nextflow"` in scripts, `Prism/bin/nextflow` by hand from
  the project folder.
  - A launcher installed with `NXF_VER` set defaults to that version
    (`NXF_VER=${NXF_VER:-'25.10.4'}`).
  - `bin/nextflow` is not tracked in git.
- **Why it is pinned:** Nextflow 26.04 made strict syntax the default,
  which broke older scripts.
  - Check the pipeline's `manifest.nextflowVersion`.
  - oncoanalyser 3.0.0 has `'!>=25.10.4'`, with no maximum.
- **Nextflow 26.10:** `executor.cpus` and `executor.memory` also cap SLURM
  jobs. Do not set them in `executor {}`.

## 5. Checking configs and scripts

From the project folder:
```bash
bash -n Prism/runOncoanalyser.sh
shellcheck -x Prism/runOncoanalyser.sh Prism/00.SETUP.sh
Prism/bin/nextflow lint Prism/conf/iris.config
Prism/bin/nextflow -c Prism/conf/iris.config config -flat
```
- **`shellcheck`** reports SC1091 for `source "$PDIR/conf/ENV.sh"`; that
  is expected.
- **Syntax:** `nextflow lint <file>`. The user said "lint"; do not read
  it as "link".
- **Resolved values:** `nextflow -c <file> config -flat`. `-c` goes
  before `config`; `nextflow config -c <file>` is wrong.
- **If `NXF_HOME` does not exist** (the user may have deleted it), do not
  recreate the user's folder. Use a throwaway one:
  `NXF_VER=25.10.4 NXF_HOME=/localscratch/soccin/<task>/nxf_home Prism/bin/nextflow lint <file>`.
  Both commands download the runtime into `NXF_HOME` the first time.
- **"Write X, do not run"** means: write it, check it as above, and stop.
  Do not execute it or submit it.
- Do not start a pipeline run or submit a job unless the user asks.

## 6. SLURM and sbatch

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
- `$0` inside the job is SLURM's copy of the script. See the next
  subsection for how to find the repo.
- Do not pass `--export=` with anything other than `ALL`. Since April 2026
  that triggers a broken user-environment probe on isca001 / bic_devs.
- Nextflow's task submissions also go through this wrapper. This works:
  `scontrol write batch_script` for task 10652310 showed the inserted
  exports.

### Finding the repo inside a job
`sbatch` runs a copy of the script from the SLURM spool folder, so `$0`
does not point into the repo there. Prism's run script tries these in
order and uses the first folder that holds a marker file
(`conf/iris.config`, checked by `is_prism_dir`):
1. `dirname "$(readlink -f "$0")"`: a direct run.
2. `$SBATCH_SCRIPT_DIR`: set by the wrapper.
3. The `Command=` field of `scontrol show job $SLURM_JOB_ID`: the
   submitted script's path, for plain `/usr/bin/sbatch`, which does not
   set `SBATCH_SCRIPT_DIR`. Only tried when `SLURM_JOB_ID` is set.

If none matches, it stops with an error. The marker check keeps a wrong
candidate (the spool folder, or a stale `SBATCH_SCRIPT_DIR`) from being
used. Steps 2 and 3 have not been tested in a job (section 14).

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
- **Choose by the time each task requests,** as in `conf/iris.config`:
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

## 7. Run script conventions (see runOncoanalyser.sh)

- **Arguments:**
  - The input samplesheet is a required first argument. No
    argument-free commands.
  - The optional second argument is the work folder of a run to resume
    (section 9).
  - Keep `input` out of `params.yaml`. A command-line `--input` overrides
    the params file anyway.
- **Help and errors:**
  - `-h` / `--help` prints usage and exits before anything else runs.
  - A wrong argument count prints usage to stderr and exits 1.
- **Checks:** `set -eu`. Check that the input file exists (`realpath`),
  and that the optional `WORK_DIR` exists.
- **Repo folder:** find `PDIR` as in section 6 and stop with an error if
  it is not found.
- **Nextflow:** stop if `$PDIR/bin/nextflow` is missing, and tell the
  user to run `$PDIR/00.SETUP.sh`.
- **It must work both via `sbatch` and run directly.** Nothing in the
  head script may depend on `SLURM_*` variables. The one use, the
  `scontrol` lookup, is a fallback that runs only when `SLURM_JOB_ID` is
  set.
- **Paths:**
  - `source "$PDIR/conf/ENV.sh"`; do not copy its exports into the
    script.
  - Every repo file by `$PDIR/...`.
  - The samplesheet as an absolute path.
  - Outputs relative to the folder the script is started from.
- **Output:** print at the start:
  - `PRISM=$PDIR (<git describe --always --dirty>)`, which records which
    wrapper commit ran
  - `INPUT`, `NXF_WORK` and `TMPDIR`
- **Header:** `-J <pipeline>` and `-o <pipeline>.%j.out`, plus the head
  job resources in section 6. The `-o` path is relative, so the log goes
  to the submit folder.
- **Finish:** `chmod +x`, then check as in section 5.
- **Usage lines:**
  - `sbatch ./Prism/runOncoanalyser.sh INPUT.csv [WORK_DIR]`
  - `./Prism/runOncoanalyser.sh INPUT.csv [WORK_DIR]`

## 8. 00.SETUP.sh

- Run once per clone, from the project folder:
  `./Prism/00.SETUP.sh`. It is executable (mode 755), as in adagio.
  - What went wrong: Prism first committed it as 644, so it could only
    be run with `bash`.
- It sources `conf/ENV.sh`. That keeps `NXF_VER` in one file, and the
  launcher's first download of the Nextflow runtime goes to `NXF_HOME`,
  not `~/.nextflow`.
- It follows the `/tmp` rule like the run script: a
  `mktemp -d -p /localscratch/soccin/<pipeline> setup_XXXXXX` folder,
  with `TMPDIR` and `-Djava.io.tmpdir` pointing at it.
- It then runs `curl -s https://get.nextflow.io | bash` in `$PDIR/bin`.

## 9. Work folders and -resume

- **A new work folder per run:** `NXF_WORK` gets a new timestamped folder
  each time `conf/ENV.sh` is sourced.
- **`-resume` only reuses tasks from the same folder,** so the run script
  takes the old folder as an optional second argument. Also:
  - Resume with the same samplesheet.
  - `-resume` needs `.nextflow/` (the run history) in the launch folder,
    so resume from the same project folder. Without it, Nextflow prints
    "Option `-resume` is ignored".
- **Large remote files** go into the work folder and can be cut off.
  - Example: the 27 GB STAR index tarball. Nextflow retries the whole
    file.
  - Download big references yourself instead, with `curl -C -` in an
    sbatch job.

## 10. Config files

Each file has one job:
- `conf/ENV.sh`: environment variables only.
- `conf/iris.config`: infrastructure only (executor, queue, limits,
  retries, Singularity, temp setup). No pipeline parameters.
- `conf/genome_<genome>.config`: reference paths.
- `conf/params.yaml`: `mode`, `genome`, `outdir`. No `input`.

**oncoanalyser `errorStrategy`:** its tools exit 1 on real failures,
hence `{ task.exitStatus != 1 ? 'retry' : 'finish' }` with
`maxRetries = 3`. Check which exit codes a new pipeline's tools use
before copying this.

Values that appear in more than one file and must agree:
- The Singularity cache path: three lines in `conf/ENV.sh` and
  `singularity.cacheDir` in `conf/iris.config`.
- The temp root `/localscratch/soccin/<pipeline>`: the run script,
  `00.SETUP.sh`, and twice in `conf/iris.config` (`beforeScript` and
  `runOptions`).
- The pipeline version: `-revision` in the run script and the reference
  data path in the genome config.

## 11. Adapting Prism to a new pipeline

Copy Prism's files into the new repo, then change:
- **Names:** the run script's file name, `#SBATCH -J` and `-o`, the
  header comment and usage text, the `PRISM=` label, and
  `is_prism_dir` (and its marker, if the config file is renamed).
- **Temp root** `/localscratch/soccin/oncoanalyser`: in the run script,
  `00.SETUP.sh`, and `beforeScript` and `runOptions` in
  `conf/iris.config`.
- **`conf/ENV.sh`:** `NXF_HOME` and `NXF_WORK` go under
  `/scratch/core001/bic/socci/<Pipeline>/`. Set `NXF_VER` from the new
  pipeline's `manifest.nextflowVersion` (section 12). Keep the three
  cache lines.
- **`conf/iris.config`:** check `errorStrategy` against the new tools'
  exit codes, and `resourceLimits` against the new pipeline's largest
  process. The queue closure, the `executor` block and the Singularity
  settings carry over.
- **Run script:** the pipeline name, `-revision`, and the config and
  params files it passes.
- **Pipeline-specific files to rewrite:** `conf/params.yaml`, the genome
  config, `scripts/`, `examples/`.
- **Docs:** `README.md`, `CLAUDE.md`, and the pipeline instructions in
  `docs/`. Copy this file and add the new pipeline's lessons to it.

Then grep the new repo for the old names, for example
`grep -rniE 'oncoanalyser|prism|wigits|grch38_hmf' <repo>`, and check
that every hit is meant to be there.

## 12. nf-core facts to check for each pipeline

From oncoanalyser 3.0.0; look up the new pipeline's values.
- **`manifest.nextflowVersion`** in `nextflow.config` at the release tag.
  A leading `!` makes the minimum a hard limit. The `NXF_VER` values in
  the pipeline's CI show which Nextflow versions are tested.
- **Label times decide the partition.** The queue closure uses
  `task.time`, which comes from the labels in the pipeline's
  `conf/base.config`. oncoanalyser: `process_single` and `process_low`
  4h, `process_medium` 8h, `process_high` 16h, `process_long` 20h, all
  multiplied by `task.attempt`. So no task goes to `cpushort` or
  `cmobic_short` unless `withName:` blocks set shorter times.
- **`check_max(...)`:** the nf-core docs' retry example uses this old
  helper. With `resourceLimits`, write `memory = { 64.GB * task.attempt }`
  instead; do not copy the `check_max` line.
- **`-c` is for infrastructure,** not pipeline parameters, which go in
  `-params-file`. The genome config, which sets `params.genomes`, is the
  documented exception.

## 13. Docs

- After every script or config change, update `README.md` and the
  pipeline instructions doc in `docs/` to match.
- **Do not paste whole scripts into docs.** Describe the steps and name
  the file.
  - What went wrong: on 2026-11-14 the oncoanalyser instructions doc
    still held a full copy of the run script from before the move to
    Prism. It no longer matched the script and was replaced by a list of
    steps.
- **After moving or renaming files, grep every doc for the old paths.**
- **Do not write claims that no file or test supports.**
  - What went wrong: the first draft of Prism's `CLAUDE.md` said a run
    takes days and uses many nodes, and that new pipelines would get
    their own run script in Prism. Neither came from a file. Both were
    removed, and the second was wrong: new pipelines get their own repo.

## 14. Status

- **2026-11-13, project-folder layout:** the four checks in section 2
  passed, on work folder
  `/scratch/core001/bic/socci/Oncoanalyser/20261113_131345`. The
  configs moved into Prism with only comment changes.
- **Not checked:** what is actually in a task node's `/localscratch`.
  That needs a shell on that node, and `srun` is not allowed.
- **Not run in the Prism layout (as of 2026-11-14):**
  - `00.SETUP.sh`
  - a direct run
  - a run through the `~/bin/sbatch` wrapper (`SBATCH_SCRIPT_DIR`)
  - a run through `/usr/bin/sbatch` (the `scontrol` lookup)

  A direct test run was blocked by the session's permission check, so
  the checks stayed at `bash -n`, `shellcheck -x` and `-h`.
- Record results here when these are run.

## 15. Working with this user

- **Files:**
  - Do not edit files you did not create (the CLAUDE.md files, the
    user's project `ENV.sh`, `~/bin/sbatch`) unless asked. Propose the
    text instead.
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
- **Spelling:** the pipeline is "oncoanalyser" (a proper name); use US
  spelling for everything else. No emoji.
- **Certainty:** ask until you are 95% sure. Do not guess paths,
  partitions or limits.

## 16. Checklist for a new pipeline

1. Read both CLAUDE.md files, this file, and adagio. Run `git log` in the
   new repo.
2. Ask the user what you cannot read from files: the repo name, where
   the run script goes (top or `bin/`), the pipeline version and mode,
   and later, public or private.
3. Copy Prism's files and change everything in section 11. Grep for the
   old names.
4. In `conf/ENV.sh`, put `NXF_HOME` and `NXF_WORK` under
   `/scratch/core001/bic/socci/<Pipeline>`, reuse the shared image
   cache, and pin `NXF_VER` from `manifest.nextflowVersion`. Nothing
   goes in `~`.
5. Add `00.SETUP.sh` (section 8) and the `.gitignore` lines (section 1).
6. Add the head temp setup to the run script and `00.SETUP.sh`, and the
   task temp setup (`beforeScript` and the two `-B` mounts) to the
   config (section 2).
7. Add the queue closure, `resourceLimits` with time `7.d`, and an
   `errorStrategy` that fits the new tools' exit codes.
8. Write the run script as in section 7, with the repo lookup from
   section 6.
9. Check with `bash -n`, `shellcheck -x`, `nextflow lint` and
   `nextflow -c ... config -flat` (section 5).
10. Write `README.md`, `CLAUDE.md` and the instructions doc (section 13).
11. Do not run anything. Tell the user what to test: `00.SETUP.sh`, a
    direct run, a run through the wrapper, a run through
    `/usr/bin/sbatch`, and the first task's `.command.run` (section 2).
12. Commit in logical groups only when asked. Push only after "Are you
    sure?".
