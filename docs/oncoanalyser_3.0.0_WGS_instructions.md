# Running nf-core/oncoanalyser 3.0.0 on WGS data

Source: https://nf-co.re/oncoanalyser/3.0.0/docs/usage/ and
https://nf-co.re/oncoanalyser/3.0.0/docs/usage/compute_resources/
(reviewed 2026-11-12). Commands below are adapted to IRIS (SLURM,
Singularity). They are run from the project folder, with Prism cloned
into it as `Prism/`. File paths without a `Prism/` prefix are relative
to the Prism folder.

## 1. What the pipeline does in WGS mode

`--mode wgts` runs the Hartwig WiGiTS toolchain on whole genome DNA
(and optionally whole transcriptome RNA):

- Alignment (BWA-MEM2, only if you start from FASTQ)
- REDUX (duplicate marking, BQR, jitter/MS tables)
- AMBER, COBALT (BAF, read depth)
- SAGE (somatic + germline SNV/indel), PAVE (annotation)
- ESVEE (structural variants)
- PURPLE (purity, ploidy, copy number)
- LINX, CHORD, LILAC (HLA), CUPPA, VIRUSBREAKEND/VIRUSINTERPRETER,
  TEAL, PEACH, ORANGE (summary report)

Supported sample setups: tumor/normal DNA, tumor-only DNA, with or
without tumor RNA, and tumor/normal/donor.

## 2. Prerequisites

### 2a. Nextflow version

- `nextflow.config` at tag 3.0.0 declares
  `nextflowVersion = '!>=25.10.4'`. The `!` makes this a hard
  minimum: older versions refuse to run the pipeline.
- No maximum version is declared.
- 3.0.0 was migrated to Nextflow strict syntax (release note
  "Migrate to strict syntax", PR #312), so the strict parser that
  became the default in Nextflow 26.04 does not break it.
- The pipeline's CI tests `NXF_VER` values `25.10.4` and
  `latest-everything`. The latest-version job uses
  `continue-on-error`, so a failure there does not fail CI.
  3.0.0 was released Sep 1, before Nextflow 26.10.0 (Oct 12), so
  "latest" at release time was 26.04.x.
- Nextflow 26.10 made no syntax changes. Its breaking changes that
  matter on SLURM: `executor.cpus` and `executor.memory` now also
  limit grid executor tasks, and tasks requesting more fail. Do not
  set those in `executor {}`. 26.10 is untested with oncoanalyser.

Recommendation: use **25.10.4**, the version the pipeline's CI
gates on. 26.04.x should also work. Avoid 26.10.x for now.

### 2b. Local Nextflow install in `Prism/bin/`

Install the pinned launcher once per clone, and always call
`Prism/bin/nextflow`:

```bash
./Prism/00.SETUP.sh
```

`00.SETUP.sh` sources `conf/ENV.sh` (which sets `NXF_VER=25.10.4`)
and runs `curl -s https://get.nextflow.io | bash` in `Prism/bin`.
The launcher installed with `NXF_VER` set keeps that version as its
default: `bin/nextflow` contains `NXF_VER=${NXF_VER:-'25.10.4'}`
(checked 2026-11-12). `bin/nextflow` is not tracked in git.

### 2c. Keep everything out of the home folder

Nextflow and Singularity write under `~` (and `work/` under the
launch directory) by default. These are set in `conf/ENV.sh`, which
`runOncoanalyser.sh` and `00.SETUP.sh` source. Source it yourself
before calling `Prism/bin/nextflow` by hand.
Pipeline folders live under `/scratch/core001/bic/socci/Oncoanalyser/`;
the Singularity cache is the existing shared
`/scratch/core001/bic/socci/singularity-cache`.

```bash
export NXF_VER=25.10.4
export NXF_HOME=/scratch/core001/bic/socci/Oncoanalyser/nextflow_home
export NXF_WORK=/scratch/core001/bic/socci/Oncoanalyser/$(date +%Y%m%d_%H%M%S)
export NXF_SINGULARITY_CACHEDIR=/scratch/core001/bic/socci/singularity-cache
export SINGULARITY_CACHEDIR=/scratch/core001/bic/socci/singularity-cache
export APPTAINER_CACHEDIR=/scratch/core001/bic/socci/singularity-cache
export NXF_OPTS='-Xms1g -Xmx4g'
```

| Variable | Default | Holds |
|---|---|---|
| `NXF_HOME` | `~/.nextflow` | Nextflow runtime jars, pulled pipeline code (`assets/`), plugins (nf-schema) |
| `NXF_WORK` | `./work` in the launch directory | Task work directories (intermediate files; several TB per WGS pair) |
| `NXF_SINGULARITY_CACHEDIR` | none (`singularity.cacheDir`) | Pipeline container images |
| `SINGULARITY_CACHEDIR` / `APPTAINER_CACHEDIR` | `~/.singularity` / `~/.apptainer` | Singularity's own download cache |
| `NXF_OPTS` | none | Limits the Nextflow head process JVM memory |

`.nextflow/` and `.nextflow.log` are written to the launch directory:
launch from the project directory, never from `~`.

### 2d. Other requirements

- A container engine. On this cluster use `-profile singularity`
  (or `apptainer`). Docker is not an option on shared HPC.
- Disk: REDUX and alignment each request up to 750 GB scratch per
  task in Hartwig's settings; plan for several TB of `work/` space
  per tumor/normal pair. `NXF_WORK` puts it on shared storage
  (`/scratch/core001/bic/socci/Oncoanalyser/<YYYYMMDD_HHMMSS>`), not
  `/localscratch` or the project directory, because SLURM tasks run
  on other nodes.

## 3. Choose the genome

Use one of the bundled Hartwig genomes:

- `GRCh38_hmf` (recommended for new work)
- `GRCh37_hmf`

Important constraints if you start from BAM/CRAM:

- The BAMs **must be aligned to the bundled genome** (same FASTA,
  same contigs).
- Supplementary reads must be soft-clipped (BWA-MEM2 `-Y`) and mate
  CIGAR tags (MC) present. Otherwise run Picard FixMateInformation.
- For GRCh38: HLA class I ALT contigs break LILAC and variant
  calling. Remove or hard-mask them, or realign with Bamtools
  AltContigRemapper.

If existing BAMs were produced by a different pipeline/reference
(e.g. a b37 or GRCh38 build with ALTs), the safest route is to start
from FASTQ and let oncoanalyser align. Check the BAM header (`@SQ`
lines and `@PG` for `-Y`) before deciding.

## 4. Stage reference data once

Without staging, every run re-downloads reference data from
`https://data.oncoanalyser.com/r2/`. Stage it once to shared storage.

Reference type values for `--ref_data_types`:

| Value | Needed for |
|---|---|
| `wgs` | WGS from alignments: fasta, fai, dict, img, hmftools, gridss_index |
| `bwamem2_index` (alias `dna_alignment`) | DNA FASTQ input |
| `wts` | WTS from alignments: fasta, fai, dict, img, hmftools |
| `star_index` (alias `rna_alignment`) | RNA FASTQ input |
| `targeted`, `panel` | Targeted/panel analysis. The only built-in panel, TSO500, is GRCh37-only. |

Minimum sets: WGS from FASTQ needs `wgs,bwamem2_index`; WGS from BAM
needs `wgs`; RNA from BAM adds `wts`; RNA from FASTQ adds `star_index`.

This project stages the full GRCh38 suite for WGS and WTS from any
input type, so the choice of FASTQ vs BAM and DNA vs RNA can be made
later without another download:

```
--ref_data_types wgs,wts,bwamem2_index,star_index
```

`targeted`/`panel` are left out: there is no built-in GRCh38 panel.

Run this as an `sbatch` job (it downloads and unpacks large files),
after `source Prism/conf/ENV.sh`. For this project the outdir is
`/scratch/core002/nds/nf-core`:

```bash
Prism/bin/nextflow run nf-core/oncoanalyser \
  -revision 3.0.0 \
  -profile singularity \
  -c Prism/conf/iris.config \
  --mode prepare_reference \
  --ref_data_types wgs,wts,bwamem2_index,star_index \
  --genome GRCh38_hmf \
  --outdir /scratch/core002/nds/nf-core
```

Incremental staging: types can be added later by running
`prepare_reference` again into the same outdir with only the missing
types. This project first staged `wgs,bwamem2_index`. That already
includes everything `wts` needs (FASTA, indexes, and the
`hmf_pipeline_resources` bundle with its `rna/` Isofox files), so the
only missing piece for RNA is the STAR index:

```bash
Prism/bin/nextflow run nf-core/oncoanalyser \
  -revision 3.0.0 \
  -profile singularity \
  -c Prism/conf/iris.config \
  --mode prepare_reference \
  --ref_data_types star_index \
  --genome GRCh38_hmf \
  --outdir /scratch/core002/nds/nf-core
```

Staged layout (`/scratch/core002/nds/nf-core/reference_data/3.0.0/`):

```
GRCh38_masked_exclusions_alts_hlas.fasta (+ .fai, .dict, .img)
bwa-mem2_index-2.2.1/
gridss_index-2.13.2/
hmf_pipeline_resources.38_v3.0.0--8/   (common/ dna/ misc/ rna/)
star_index-gencode_38-2.7.3a/
```

The STAR tarball is 27.4 GB and the pipeline's download can be cut
off ("Read data length does not match expected size"). If all 3
retries fail, download it with a resumable tool (`curl -C -`) into
the folder above and `tar -xzf` it.

Output lands in `<outdir>/reference_data/3.0.0/`. Delete that run's
task folders under `NXF_WORK` afterwards. Alternative: download the
tarballs listed on the usage page manually and `tar -xzvf` them.

`conf/genome_GRCh38_hmf.config` points the pipeline at the staged
files (passes `nextflow lint`). It sets `params.genomes.GRCh38_hmf`
(`fasta`, `fai`, `dict`, `img`, `bwamem2_index`, `gridss_index`,
`star_index`) and `ref_data_hmf_data_path` to the paths above. Use it
with `--genome GRCh38_hmf`.

## 5. Build the samplesheet

CSV with header. Columns:

| Column | Meaning |
|---|---|
| `group_id` | One analysis unit; top-level output dir |
| `subject_id` | Same for all rows in a group |
| `sample_id` | Tumor ID prefixes most output files |
| `sample_type` | `tumor`, `normal`, `donor` |
| `sequence_type` | `dna`, `rna` |
| `filetype` | `fastq`, `bam`, `bai`, `cram`, `crai`, `redux_dir`, ... |
| `info` | `key:value;key:value` (FASTQ needs `library_id` and `lane`) |
| `filepath` | Absolute path; for FASTQ, `R1;R2` |

Optional `info` keys: `cancer_type` (used by NEO and the ORANGE
header), `read_group_overrides` (e.g. `ID=id-1|SM=sample-1|LB=lib-1`).

### 5a. FASTQ, tumor/normal (one row per lane per sample)

Only gzipped, non-interleaved paired-end FASTQ.

```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath
P1,P1,P1-N,normal,dna,fastq,library_id:N1;lane:001,/abs/P1-N_L001_R1.fastq.gz;/abs/P1-N_L001_R2.fastq.gz
P1,P1,P1-N,normal,dna,fastq,library_id:N1;lane:002,/abs/P1-N_L002_R1.fastq.gz;/abs/P1-N_L002_R2.fastq.gz
P1,P1,P1-T,tumor,dna,fastq,library_id:T1;lane:001,/abs/P1-T_L001_R1.fastq.gz;/abs/P1-T_L001_R2.fastq.gz
P1,P1,P1-T,tumor,dna,fastq,library_id:T1;lane:002,/abs/P1-T_L002_R1.fastq.gz;/abs/P1-T_L002_R2.fastq.gz
```

(The docs show only a tumor-only FASTQ example; the normal rows above
follow the same format.)

### 5b. BAM, tumor/normal

```csv
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,filepath
P1,P1,P1-N,normal,dna,bam,/abs/P1-N.dna.bam
P1,P1,P1-T,tumor,dna,bam,/abs/P1-T.dna.bam
```

`.bai` must sit next to the BAM as `<bam>.bai`; otherwise add a `bai`
row. Tumor-only: drop the normal row. Multiple patients: one
`group_id`/`subject_id` per patient in the same sheet.

### 5c. CRAM

Use `cram`/`crai`. CRAM must be version <= 3.0, with NM tags stored.
Small CRAMs (< ~10 GB) read >10x slower than BAM; converting to BAM
first is faster.

## 6. Cluster config (`conf/iris.config`)

The cluster settings live in `conf/iris.config`. `runOncoanalyser.sh`
passes it with `-c`. Per the docs, `-c` is for resources and
executor settings, not pipeline parameters (the reference config in
step 4 is the documented exception).

What it sets:

- `executor = 'slurm'`
- Partition chosen from the time each task requests (`task.time`),
  re-evaluated on each retry:

  | Requested time | Partition | Partition walltime limit |
  |---|---|---|
  | <= 2h | `cpushort,cmobic_short` | cpushort 2h |
  | > 2h and <= 3h | `cmobic_short` | 3h |
  | > 3h, or none | `cmobic_cpu` | 7 days |

- QOS `priority` (`--qos=priority` through `clusterOptions`) for
  tasks on `cmobic_cpu` and `cmobic_short`, none for tasks listed for
  `cpushort`, which allows only `normal`. `priority` puts the tasks
  ahead of `normal` jobs waiting for the cmobic partitions; it has no
  limits and does not change preemption.
- `resourceLimits = [ cpus: 52, memory: 975.GB, time: 7.d ]`: the most
  SLURM can allocate on any node of the three partitions (4 cores and
  20 GB per node are reserved, and `cpushort` sets
  `MaxCPUsPerNode=52` and `MaxMemPerNode=998557` MB). The time cap is
  the `cmobic_cpu` limit.
- `withName: 'PURPLE' { memory = { 96.GB * task.attempt } }`. PURPLE
  loads the whole somatic VCF into the Java heap, and the default
  36 GB failed with `OutOfMemoryError` on a tumor-only sample.
- ESVEE: 24 CPUs (not scaled by attempt) and `216.GB * task.attempt`
  (default 12 CPUs / 72 GB). See the notes below.
- SAGE somatic: `24 * task.attempt` CPUs and `144.GB * task.attempt`
  (default 12 CPUs / 72 GB), the pipeline's own attempt-2 size.
- Retry every failure once; after that, retry any exit code except 1,
  up to 3 retries in all.
- `executor { queueSize = 200; submitRateLimit = '10 sec' }`. At 50,
  alignment used about half of `cmobic_cpu`, which fits about 91
  12-CPU tasks; 200 fills it, and the rest wait in SLURM with QOS
  `priority`. Not yet tested at 200 (see `ISSUES.md`).
- Singularity enabled, `autoMounts`, `cacheDir` on shared storage,
  `pullTimeout = '2h'`.
- No writes to the node's `/tmp`. Each task (one SLURM job) gets its
  own folder `/localscratch/soccin/oncoanalyser/<jobid>` on the node.
  `beforeScript` creates it and points `TMPDIR`, `SINGULARITY_TMPDIR`
  and `APPTAINER_TMPDIR` at it, replacing the head job's values that
  sbatch passes to every task. `singularity.runOptions` mounts the
  folder twice: over the container's `/tmp` (for tools that write to
  `/tmp` directly) and at its own path (so `TMPDIR` exists inside the
  container). There is no automatic cleanup; remove the folders by hand.
- Before a full run, check the `singularity exec` line in the first
  task's `.command.run` contains
  `-B /localscratch/soccin/oncoanalyser/$SLURM_JOB_ID:/tmp -B /localscratch/soccin/oncoanalyser/$SLURM_JOB_ID`
  and the `beforeScript` block sets `TMPDIR`.

`cacheDir` must match `NXF_SINGULARITY_CACHEDIR` in `conf/ENV.sh`
(section 2c).
Do not add `cpus` or `memory` to the `executor {}` block: since
Nextflow 26.10 those cap SLURM tasks too.

Notes:

- With oncoanalyser's default labels (`conf/base.config`) every task
  requests at least 4h (`process_single`/`process_low` 4h,
  `process_medium` 8h, `process_high` 16h, `process_long` 20h, all
  multiplied by `task.attempt`), so nothing goes to the short
  partitions until short tools get `withName:` time overrides. Set
  those from real runtimes in the trace report of a first run.
- Set `resourceLimits` cpus and memory to what SLURM can allocate on
  the smallest node across the partitions a task can be sent to, so
  every request fits on any of them. A job listed for several
  partitions must be valid in each one (`EnforcePartLimits=ALL`).
  Check a request with `/usr/bin/sbatch --test-only`, which submits
  nothing.
- WiGiTS tools exit with code 1 on any real failure, and a Java
  `OutOfMemoryError` also exits 1. So every failure is retried once,
  with memory doubled by `task.attempt`; after that only non-1 exit
  codes (SLURM OOM kill, node issues) are retried. A real tool failure
  therefore runs twice before the pipeline stops. Before this rule, a
  PURPLE and an ESVEE heap error each stopped a whole run.
- ESVEE memory depends on the number of chimeric reads, not the
  depth. In a run of 15 WGS groups (Nov 2026), on the default 72 GB,
  every group with more than about 300 M chimeric reads (`pctChimeras`
  x total reads, tumor plus normal) was killed at the prep BAM sort
  (exit 137) or ran out of Java heap (exit 1). One with 571 M chimeric
  reads was still killed with 144 GB and passed with 216 GB. Groups
  with few chimeric reads finished ESVEE in 10 to 75 min, so the larger
  request costs little for them. 24 CPUs ran ESVEE prep 1.4 to 2 times
  faster than 12; 36 was no faster than 24. A resume with the current
  settings (2026-11-18) ran the other large group (601 M chimeric
  reads, 918 M prep BAM records) on 24 CPUs / 216 GB: passed on
  attempt 1 in 3h56m, with no slowdown from a full heap. A group
  larger than these may still be killed at the sort on 216 GB; it is
  then retried on 432 GB.
- `MaxRSS` (sacct) equals the request for every task that reads and
  writes a lot of data, Java or not: also for alignment (bwa-mem2 and
  sambamba, trace `peak_rss` 24 to 29 GB of 72), AMBER, COBALT and
  BAMTOOLS. Do not size from it. The trace's `peak_rss` is usable,
  except for Java tools whose heap fills to the request (SAGE) and for
  ESVEE, where it is above the request.
- SAGE somatic used 97% of 12 CPUs for 1 to 5.3 h. The gain from
  24 CPUs is expected but not yet measured: every SAGE somatic task
  was cached in the 2026-11-18 resume. See `ISSUES.md`.
- ESVEE asks for 32 CPUs / 96 GB in Hartwig's settings. That is less
  memory than the 144 GB on which the group above was still killed.
- The docs' retry example uses `check_max(...)`, an old nf-core
  helper. With `resourceLimits` write `memory = { 64.GB * task.attempt }`
  instead; do not copy the `check_max` line.
- Hartwig's per-process recommendations (from the compute resources
  page) can be added as `withName:` blocks if defaults are too small:
  REDUX 72 GB/16 CPU, SAGE 72 GB/16 CPU, ESVEE 96 GB/32 CPU,
  LILAC 64 GB/16 CPU, PURPLE 40 GB/8 CPU, alignment 72 GB/12 CPU.
  Use the selector from the pipeline's `conf/modules.config` (for
  example `'.*:ESVEE_CALLING:ESVEE'`), which matches the full process
  name `NFCORE_ONCOANALYSER:WGTS:ESVEE_CALLING:ESVEE`.
- Pre-pulling the Singularity images into `cacheDir` once avoids
  timeouts on the first run.

## 7. Run the analysis

The fixed run parameters are in `conf/params.yaml`. The samplesheet is
not in it; it is passed to the script as a required argument.

```yaml
mode: 'wgts'
genome: 'GRCh38_hmf'
outdir: './out'
```

The Nextflow head process runs from `runOncoanalyser.sh`, either as its
own SLURM job (long-lived, low resource) or in the current shell. Start
it from the project directory. Do not run a full run on isca006 and do
not use `srun`. `./Prism/runOncoanalyser.sh -h` prints the usage.

New run (new work folder):

```bash
sbatch ./Prism/runOncoanalyser.sh oncoanalyser_input.csv
```

Resume a previous run (reuse its work folder):

```bash
sbatch ./Prism/runOncoanalyser.sh oncoanalyser_input.csv /scratch/core001/bic/socci/Oncoanalyser/<YYYYMMDD_HHMMSS>
```

The script converts INPUT.csv to an absolute path and passes it as
`--input`. Use the same samplesheet when resuming; a different one
changes the task inputs and those tasks rerun.

On resume, retries start again at attempt 1. If a task failed for
lack of memory, raise its starting memory in `conf/iris.config`
first. Resource, partition, QOS and `errorStrategy` settings are not
part of the task hash, so changing them keeps finished tasks cached.
`ext.*` settings used by a module script are part of it: changing one
reruns that process and everything after it. In the first minutes,
`.nextflow.log` should show the finished tasks as `Cached process`.

What `runOncoanalyser.sh` does (read the script for details):

1. Checks the arguments and the samplesheet.
2. Finds the Prism folder. `sbatch` runs a copy of the script from the
   SLURM spool folder, so it tries the folder of the script,
   `SBATCH_SCRIPT_DIR` (set by the `~/bin/sbatch` wrapper), and the
   submitted path from `scontrol show job`, in that order.
3. Checks that `Prism/bin/nextflow` is installed.
4. Sources `conf/ENV.sh` and, if given, sets `NXF_WORK` to `WORK_DIR`.
5. Sets up the head job temp folder (below).
6. Prints `PRISM=` (Prism folder and git version), `INPUT=`,
   `NXF_WORK=` and `TMPDIR=`.
7. Runs:

   ```bash
   "$PDIR/bin/nextflow" run nf-core/oncoanalyser \
       -revision 3.0.0 \
       -profile singularity \
       -config "$PDIR/conf/genome_GRCh38_hmf.config" \
       -c "$PDIR/conf/iris.config" \
       -params-file "$PDIR/conf/params.yaml" \
       --input "$INPUT" \
       -resume
   ```

`conf/ENV.sh` sets `NXF_WORK`, so no `-work-dir` is needed; `work/` goes
to `/scratch/core001/bic/socci/Oncoanalyser/<YYYYMMDD_HHMMSS>`, a new
folder each time `conf/ENV.sh` is sourced (e.g. `20261112_193000`).
`-resume` only reuses tasks from the same work folder, so to resume
pass that folder as the second argument. The samplesheet and folder in
use are printed in `oncoanalyser.<jobid>.out` as `INPUT=...` and
`NXF_WORK=...`; the folder is also shown as `workDir` in the run header
and in `.nextflow.log`. Results go to `./out` in the project directory.

The head job writes its temporary files (Nextflow JVM, Singularity image
pulls) to a new `/localscratch/soccin/oncoanalyser/head_XXXXXX` folder (made by
`mktemp`, so it works with `sbatch` or run directly) on its node, never
`/tmp`; the path is printed as `TMPDIR=...` in the job log. No automatic
cleanup.

`-t 7-00:00:00` is the `cmobic_cpu` walltime limit; the head job
must outlive the whole run. If a run needs longer, resubmit with the
work folder as the second argument so `-resume` reuses completed tasks;
resubmitting without it starts a new folder and reruns everything.

To collect the reports for delivery, run from the project directory:

```bash
./Prism/deliverLocal.sh [OUT_DIR] [RESULTS_DIR]
```

It hard links each group's ORANGE PDF (`orange/*.orange.pdf`), QSEE
QC report (`qsee/*.qsee.vis.report.pdf` and `.png`), LINX HTML report
(`linx/*_linx.html`), PURPLE circos plot (`purple/plot/*.circos.png`)
and CUPPA plot (`cuppa/*.cuppa.vis.png`, tumor/normal groups only)
from `OUT_DIR` (default `out`) into `RESULTS_DIR/PDF`,
`RESULTS_DIR/HTML` and `RESULTS_DIR/PNG` (default `results`). Hard
links take no extra space but need both folders on the same file
system. Rerunning skips links already in place and reports, without
replacing, a different file under an existing name.

Deliverables made after the pipeline go under `res/prism/` in the
folder that holds `OUT_DIR`; their logs go to `OUT_DIR/logs/other/`,
next to the pipeline's own. To run QSEE (Hartwig's QC tool) once over
all groups:

```bash
sbatch ./Prism/qcCohort.sh [OUT_DIR]
```

oncoanalyser runs QSEE once per group (`modules/local/qsee/main.nf`
passes one `-tumor` and `-reference`, never `-sample_id_file` or
`-merge_plots`). `qcCohort.sh` copies each sample's QSEE input files
from `OUT_DIR` to a temp folder under `/localscratch`, then runs the
QSEE 1.0 image from the pipeline's Singularity cache with
`-sample_id_file` and `-merge_plots`. Reference files come from
`ref_data_hmf_data_path` in the newest
`OUT_DIR/pipeline_info/params_*.json`. Output in `res/prism/qc/hmf/`:

- `multisample.qsee.vis.report.<set>.pdf`: one page per tumor and its
  normal, against the Hartwig cohort percentiles (the same plots as
  the per-group reports; QSEE 1.0 does not plot samples together)
- `sample_ids.<set>.tsv`: the sample list given to QSEE
- `data/multisample.qsee.status.<set>.tsv.gz`: PASS/WARN/FAIL per
  sample and metric
- `data/multisample.qsee.vis.data.<set>.tsv.gz`: metric values per
  sample

QSEE writes everything to one folder and its plot step reads the
metric table back from there, so the script points it at `data/` and
then moves the PDF up. Logs in `OUT_DIR/logs/other/`:
`qsee_multisample.<set>.log`, and `qcCohort.<jobid>.out`, which the
job moves there from the submit folder at its end (sbatch only; only a
file with exactly that name is moved).

`<set>` is `tumor_normal` or `tumor_only`. QSEE pairs tumors and
normals by position in the sample list and stops unless every tumor
has a normal or none does, so the two kinds of group run separately.
A set with only one sample gets QSEE's single-sample file names. The
job asks for 4 CPUs, 16 GB and 1 hour (`cpushort,cmobic_short`); one
group took about 20 seconds in the pipeline.

To write the cohort report:

```bash
sbatch ./Prism/cohortReport.sh [-m MANIFEST] [-g GENES] [OUT_DIR]
```

`cohortReport.sh` runs `scripts/cohortReport.R` (`Rscript --vanilla`;
packages `tidyverse`, `fs`, `openxlsx`), which reads the PURPLE, LINX,
Virus Interpreter, CHORD, SIGS, LILAC, PEACH and QSEE tables of every
group and writes `res/prism/cohort/prism_report.xlsx` and one CSV per
table in `data/`. The report counts and lists events; it applies no
significance test and leaves no sample out. The first sheet is for
the genes of interest: `genes_of_interest_events` lists every event in
them, one row per event. The sheets a reader can use directly come
next (`samples`, `cohort_summary`, then the event lists), then the
supporting tables and `summary_full`. There is no per-gene count
table over the cohort; see `ISSUES.md`.
The last two sheets are `manifest`, which says what each sheet holds,
and `glossary`, which explains every column once, with the sheets
that have it. The texts of both are in `scripts/cohortReport.R`
(`manifest_text`, `glossary_text`); a column without a glossary entry
gives a warning in the log. Both are built on every run from the
sheets that are written, and they cover the sheets of their own
workbook only: `prism_sv_same_breakpoint.xlsx` has its own `glossary`
sheet. When a sheet or column is added or changed, update its text
there. Column widths are set from the content (`auto_widths`: the
95th percentile of the cell lengths, at most 50 characters); the
sheets and columns that hold long text have fixed widths
(`TEXT_WIDTHS`, `COL_WIDTHS`). The main rules:

- `cohort_summary` gives each count as `n of N (pct%)`, for all
  samples and per class. `summary_full` has the same counts as
  numbers, the groups with and without a matched normal, and the
  median, minimum and maximum of the numeric measures.
- `MANIFEST` (CSV: `sample_id`, `sample_class`, any other columns)
  sets the class of each tumor sample; without it every sample is
  `Tumor`. `cohort_summary` gives its counts for all samples and per
  class.
- `GENES` (one symbol per line, `#` starts a comment) sets the genes
  of interest. They get the first sheet.
- Without `-g` the genes of interest are the genes with a pipeline
  driver call (LINX driver catalog) of likelihood 0.8 or higher in at
  least one sample, those called in the most samples first
  (`DEFAULT_GOI_MIN_LIKELIHOOD` in the R script). The list is written
  to `res/prism/cohort/genes_of_interest_default.txt`; edit a copy and
  give it back with `-g`. A rerun without `-g` replaces the file.
- `genes_of_interest_events` writes every kind of event into the same
  columns: `event_type`, `event` (what it is), `location`, `support`
  (the evidence, as text), `caution` (weak support, shared
  breakpoint), `hmf_driver` and `driver_likelihood`. `hmf_driver` is
  the pipeline's "reported" flag under a plainer name, in every
  sheet that has it: TRUE if the event is on the pipeline's list of
  candidate driver events, which follows the Hartwig driver gene
  panel (per gene, the kinds of event that count). FALSE does not
  mean the event was not detected, and TRUE does not mean it is well
  supported or a likely driver. For copy number events the report
  sets it from the driver list (`drivers` sheet). The kinds
  are mutation, amp, partial_amp, deep_del, loh, disruption and
  fusion. A partial amplification (maximum copy number of the gene
  above 3 x ploidy, minimum not) is listed in this sheet only;
  `gene_cn` has no row for it. A fusion gives one row for each of
  its genes that is a gene of interest.
- Somatic variants are the PASS records of `purple.somatic.vcf.gz`
  with a missense, inframe, nonsense, frameshift or splice effect on
  the canonical transcript (PAVE `IMPACT` field).
- Copy number events use PURPLE's driver rules for every gene:
  amplification above 3 x sample ploidy, deep deletion below 0.5
  copies, LOH below 0.5 minor allele copies.
- Disruptions are LINX's disruptive breakends on canonical
  transcripts, with breakpoint, SV type, junction copy number and
  read support (from `purple.sv.vcf.gz`). `weak_support` marks SVs in
  a LINX cluster of category `ARTIFACT`.
- `same_bp_flag` marks a breakpoint (chromosome, position,
  orientation) found in more than one sample. They are listed in a
  workbook of their own, `prism_sv_same_breakpoint.xlsx`, not in the
  report.
- Events on chrX and chrY are removed for every sample that does not
  have both a value in a manifest column named `sex` and a matched
  normal. The sex PURPLE infers is not reported. The rule keeps the
  many X/Y events of such samples out of the tables. An event in a
  gene of the list given with `-g` is kept for every sample, in every
  sheet (`gene_cols` of `redact_xy` in the R script). Without `-g` there is
  no exception: the default list is made after the rule is applied.
- DPYD status is from PEACH, which needs the normal; samples without
  one show `not assessed`.
- The driver gene panel (for `panel_role` and the LOH rows of `gene_cn`) is read
  from `ref_data_hmf_data_path` in the newest
  `OUT_DIR/pipeline_info/params_*.json`.

The log `cohortReport.<jobid>.out` is moved to `OUT_DIR/logs/other/`
at the end of the job. A rerun replaces the workbooks and the CSVs. The
job asks for 2 CPUs, 16 GB and 1 hour; 14 groups took under a minute.

To build the cohort tables (download the region tracks once per site
first, on the login node):

```bash
./Prism/scripts/fetchTracks.sh [TRACK_DIR]
sbatch ./Prism/postOA.sh [-t TRACK_DIR] [OUT_DIR]
```

`TRACK_DIR` defaults to
`/scratch/core001/bic/socci/reference/prism_tracks/GRCh38`.
`fetchTracks.sh` downloads GIAB genome stratifications v3.6 for GRCh38
(the unions of all difficult regions and of low mappability and
segmental duplications; low mappability; segmental duplications; tandem
repeats and homopolymers; satellites; GC below 25% or above 65%; gaps;
other difficult regions and each of its parts: false and collapsed
duplications, population CNV false-positive regions, L1H, VDJ, MHC,
KIR), the ENCODE blacklist v2, and UCSC RepeatMasker and centromeres.
Each goes into its own folder; `TRACK_DIR/tracks.tsv` lists the source
URL, download date and md5 of each, and the check of each GIAB file
against the GIAB md5 list.

`postOA.sh` runs four steps: `scripts/postOA.R` (stacked LINX tables),
`scripts/svCohort.R` (the SV event table as TSV),
`scripts/tsvToParquet.py` (its typed Parquet copy) and
`docs/scripts/checkSvCohort.R` (rebuilds the sources from the table and
checks the Parquet copy). It needs `Rscript --vanilla` with `tidyverse`,
`fs`, `nanoparquet`, `GenomicRanges`, `IRanges` and `S4Vectors`, and
`python3` with `pyarrow`; it checks for `nanoparquet` and `pyarrow`
before it starts. All steps follow one rule: do not infer biology;
faithfully preserve what is there. Every assumption is checked on every
row, and a failed check stops the job.

`scripts/postOA.R`: for each kind of table it stacks the files of all
groups into one `.tsv.gz` in `res/prism/db/data/`, with `group_id`
and `sample_id` in front. The lines are copied as the tool wrote them:
no row is left out and no value is changed or reformatted. So far it
covers the LINX tables:

- somatic, every group: `svs`, `clusters`, `links`, `breakend`,
  `fusion`, `drivers`, `driver_catalog`
  (`linx/somatic_annotations/<sample_id>.linx.<table>.tsv`)
- germline, groups with a matched normal: `svs`, `breakend`,
  `disruption`, `clusters`, `links`, `driver_catalog`
  (`linx/germline_annotations/<sample_id>.linx.germline.<table>.tsv`)

LINX's `neoepitope` (input to NEO) and `vis_*` (input to the LINX
plots) files are not aggregated. Which SVs LINX reads from PURPLE, and
how its tables link to each other, is in `docs/ACTUAL_METHODS.md`
(section 4, LINX) and in the descriptions in `manifest.tsv`. The script
stops if a file is missing; it stops, and writes no further table, if
the files of a table have different headers or a line has a different
number of fields than its header. At the top of `res/prism/db/`:

- `manifest.tsv`: each file, its tool and mode, the number of groups
  and rows, the source path and what it holds
- `samples.tsv`: per group, `sample_id` (the tumor, named by the PURPLE
  purity file), `normal_id` (the other sample folder in `bamtools/`),
  `run_mode` (PURPLE) and the LINX version
- `rows_per_sample.tsv`: the rows each group adds to each table

`scripts/svCohort.R` writes `res/prism/db/sv_cohort_filtered.tsv.gz`,
the core analysis table, one row per somatic SV event; step 3 writes
the same table as `sv_cohort_filtered.parquet`, the file to analyze:

- Rows: every record of each group's PURPLE `purple.sv.vcf.gz` (FILTER
  PASS, PON and INFERRED; nothing is left out). The two breakends of an
  SV (linked by `MATEID`) form one row, with the first in the VCF as the
  start; a single or inferred breakend forms one row. The script stops
  unless the events account for every VCF record.
- Per event, without a prefix: the VCF, FORMAT and ESVEE fields that
  have the same value on both breakends of every SV (`qual`, `filter`,
  `info_SVTYPE`, `tumor_VF` and others; 36 listed in the script, 35 of
  them present in the 2-group test), stored once (checked).
- Per breakend (`start_*`, `end_*`): the other VCF fields, INFO keys
  (`info_*`) and FORMAT values of tumor and normal (`tumor_*`,
  `normal_*`); the columns of ESVEE's
  `esvee/<sample_id>.esvee.breakend.tsv` that no VCF field repeats
  (`esvee_bnd_*`, joined on VCF ID = ESVEE `Id`; inferred breakends
  have none); LINX's `breakend.tsv` rows for that end, one list element
  per transcript (`linx_bnd_*`); each region track
  (`region_*`: in the track or not, RepeatMasker name, class and family,
  ENCODE blacklist label, distance to the centromere model, Hartwig
  `sv_prep_blacklist` and `unmap_regions`); the number and list of
  other groups with a breakend at the same position, or within 50 bases,
  with the same orientation, in their PURPLE somatic, PURPLE germline
  and ESVEE unfiltered calls (`xs_*`), and the number of such
  breakends in the group's own germline file (`own_germline_*`).
- Per event: the LINX `svs.tsv` row (`linx_*`, joined on `vcfIdStart`),
  its cluster (`linx_cluster_*`), fusions, links and cluster drivers
  (lists), whether the same two-breakend event is in other
  groups (`xs_pair_*`), and derived columns (`n_breakends`, `in_linx`,
  `same_chrom`, `length`, coordinates). For some SVs whose two
  breakends share a position, LINX's start is the VCF end
  (`linx_start_is_vcf_start` FALSE); LINX's columns are kept as LINX
  wrote them.

- Unpacked for analysis, each step checked on every row: fixed pairs
  split into two columns (`AD` ref/alt, `CIPOS` and `IHOMPOS`
  lower/upper, `TPR` five/three prime span); PURPLE's `PURPLE_AF`,
  `PURPLE_CN` and `PURPLE_CN_CHANGE` reduced to each breakend's own
  value (PURPLE writes both breakends' values on each record); lists
  that always hold one value kept as one value; columns that repeat
  another column dropped (most ESVEE columns, `MATEID`). A VCF `.` is
  stored as empty.
- Lists of 1 to N values (LINX transcripts, fusions, links, drivers,
  RepeatMasker labels, group lists, VCF lists such as `ASMLNKS`, LINX
  `svs.tsv` lists such as `clusterReason`) are arrays in the Parquet
  file and joined with `|`, `,` or `;` in the TSV.

`res/prism/db/sv_cohort_filtered.columns.tsv` gives every column its
description (from the VCF header where it has one), type, list
separator and count column, and its number of non-empty values.
`scripts/tsvToParquet.py` converts every cell to that type (32-bit
integers, so R reads them as integer) and stops on any value that does
not fit. The last step, `docs/scripts/checkSvCohort.R`, undoes the
unpacking, rebuilds every PURPLE VCF record and ESVEE row, compares
them with the source files, and checks every Parquet cell against the
TSV; the job fails if any check fails. The table, its Parquet copy and
the dictionary are added to `manifest.tsv`. The exact definitions
(event pairing, joins, unpacking rules and their sources, types, lists,
coordinates, the 20 tracks, the cross-sample sets, the checks, how to
read the Parquet file in R and Python) are in `docs/prism_methods.md`.

The log `postOA.<jobid>.out` is moved to `OUT_DIR/logs/other/` at the
end of the job. A rerun replaces the files. The job asks for 2 CPUs,
32 GB and 1 hour.

## 8. Optional controls

- Skip tools: `--processes_exclude virusinterpreter,orange`
- Run only some tools:
  `--processes_manual alignment,redux,amber,cobalt,esvee,sage,pave,purple`
  You must list every upstream dependency yourself; the pipeline does
  not error if you omit one, it silently skips the downstream step.
- Restart from intermediate outputs: add `redux_dir`, `sage_dir`,
  `esvee_dir`, `purple_dir`, `amber_dir` rows to the samplesheet to
  skip completed steps without the old `work/` directory.
- UMIs (if your library prep used them): `--umi_type twist|kapa|msk|tso500`.
- Non-Illumina data: `--sequencing_platform sbx|ultima`.

## 9. Suggested order of work

1. Clone Prism into the project folder, run `./Prism/00.SETUP.sh`
   (section 2b), check the environment variables in `conf/ENV.sh`
   (section 2c), and check Singularity or Apptainer is available on
   the compute nodes.
2. Inspect the WGS inputs: FASTQ or BAM? If BAM, check reference
   contigs, `-Y`, MC tags, and HLA ALT contigs.
3. Decide GRCh37_hmf vs GRCh38_hmf.
4. Stage reference data (section 4) and check
   `conf/genome_GRCh38_hmf.config`.
5. Check `conf/iris.config` and `conf/params.yaml`, and build the
   samplesheet with `Prism/scripts/makeOncoanalyserInput.R` (see
   `README.md`).
6. Run one tumor/normal pair first; check the ORANGE report and
   PURPLE purity/ploidy before submitting the full cohort.
7. Delete the run's work folder under
   `/scratch/core001/bic/socci/Oncoanalyser/` once outputs are
   verified.

## 10. Open questions for this project

- Input format of the WGS data (FASTQ vs BAM/CRAM), and the reference
  used if BAM.
- Tumor/normal pairs or tumor-only?
- Is there RNA-seq to include (would switch to WGTS with RNA rows)?
- Target genome build for downstream comparison with other results.
- UMIs in the library prep?
