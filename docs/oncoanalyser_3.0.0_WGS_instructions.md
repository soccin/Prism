# Running nf-core/oncoanalyser 3.0.0 on WGS data

Source: https://nf-co.re/oncoanalyser/3.0.0/docs/usage/ and
https://nf-co.re/oncoanalyser/3.0.0/docs/usage/compute_resources/
(reviewed 2026-11-12). Commands below are adapted to the local SLURM
cluster. Nothing has been run.

> Copied into Prism on 2026-11-14 from the project folder where the
> pipeline was first set up. The text below uses that folder's paths,
> and the copy of `runOncoanalyser.sh` in section 7 is the version from
> that folder, not the current one. In Prism the files are:
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

### 2b. Local Nextflow install in `bin/`

Install a pinned launcher in this project's `bin/` and always call
`./bin/nextflow`:

```bash
mkdir -p bin
cd bin
export NXF_VER=25.10.4
curl -s https://get.nextflow.io | bash
cd ..
```

The launcher installed with `NXF_VER` set keeps that version as its
default: `bin/nextflow` contains `NXF_VER=${NXF_VER:-'25.10.4'}`
(checked 2026-11-12). The run scripts below still export
`NXF_VER=25.10.4`; that is harmless and makes the version explicit.

### 2c. Keep everything out of the home folder

Nextflow and Singularity write under `~` (and `work/` under the
launch directory) by default. These are kept in the project's
`ENV.sh`; source it in every script before calling `./bin/nextflow`.
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
after `source ENV.sh`. For this project the outdir is
`/scratch/core002/nds/nf-core`:

```bash
./bin/nextflow run nf-core/oncoanalyser \
  -revision 3.0.0 \
  -profile singularity \
  -c config/iris.config \
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
./bin/nextflow run nf-core/oncoanalyser \
  -revision 3.0.0 \
  -profile singularity \
  -c config/iris.config \
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

`config/genome_GRCh38_hmf.config` points the pipeline at the staged
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

## 6. Cluster config (`config/iris.config`)

The cluster settings live in `config/iris.config`. Pass it with
`-c config/iris.config`. Per the docs, `-c` is for resources and
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

- `resourceLimits = [ cpus: 32, memory: 120.GB, time: 7.d ]`; the
  time cap is the `cmobic_cpu` limit.
- Retry on any exit code except 1, up to 3 times.
- `executor { queueSize = 50; submitRateLimit = '10 sec' }`
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

`cacheDir` must match `NXF_SINGULARITY_CACHEDIR` from section 2c.
Do not add `cpus` or `memory` to the `executor {}` block: since
Nextflow 26.10 those cap SLURM tasks too.

Notes:

- With oncoanalyser's default labels (`conf/base.config`) every task
  requests at least 4h (`process_single`/`process_low` 4h,
  `process_medium` 8h, `process_high` 16h, `process_long` 20h, all
  multiplied by `task.attempt`), so nothing goes to the short
  partitions until short tools get `withName:` time overrides. Set
  those from real runtimes in the trace report of a first run.
- Set `resourceLimits` cpus and memory to the largest node in the
  partitions. ESVEE asks for 32 CPUs / 96 GB in Hartwig's settings.
- WiGiTS tools exit with code 1 on any real failure, so only non-1
  exit codes (OOM, node issues) are retried.
- The docs' retry example uses `check_max(...)`, an old nf-core
  helper. With `resourceLimits` write `memory = { 64.GB * task.attempt }`
  instead; do not copy the `check_max` line.
- Hartwig's per-process recommendations (from the compute resources
  page) can be added as `withName:` blocks if defaults are too small:
  REDUX 72 GB/16 CPU, SAGE 72 GB/16 CPU, ESVEE 96 GB/32 CPU,
  LILAC 64 GB/16 CPU, PURPLE 40 GB/8 CPU, alignment 72 GB/12 CPU.
- Pre-pulling the Singularity images into `cacheDir` once avoids
  timeouts on the first run.

## 7. Run the analysis

The fixed run parameters are in `params.yaml`. The samplesheet is not
in it; it is passed to the script as a required argument.

```yaml
mode: 'wgts'
genome: 'GRCh38_hmf'
outdir: './out'
```

The Nextflow head process runs as its own SLURM job (long-lived, low
resource) from `runOncoanalyser.sh`. Submit from the project directory.
Do not run it on isca006 and do not use `srun`. `./runOncoanalyser.sh -h`
prints the usage.

New run (new work folder):

```bash
sbatch runOncoanalyser.sh oncoanalyser_input.csv
```

Resume a previous run (reuse its work folder):

```bash
sbatch runOncoanalyser.sh oncoanalyser_input.csv /scratch/core001/bic/socci/Oncoanalyser/<YYYYMMDD_HHMMSS>
```

The script converts INPUT.csv to an absolute path and passes it as
`--input`. Use the same samplesheet when resuming; a different one
changes the task inputs and those tasks rerun.

`runOncoanalyser.sh`:

```bash
#!/bin/bash
#SBATCH -J oncoanalyser
#SBATCH -p cmobic_cpu
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 7-00:00:00
#SBATCH -o oncoanalyser.%j.out

# Nextflow head job for nf-core/oncoanalyser 3.0.0.
# Submit from the project directory:
#   sbatch runOncoanalyser.sh INPUT.csv [WORK_DIR]
# INPUT.csv: the oncoanalyser samplesheet (e.g. oncoanalyser_input.csv).
# Optional WORK_DIR: an existing
#   /scratch/core001/bic/socci/Oncoanalyser/<YYYYMMDD_HHMMSS>
# folder to -resume a previous run; otherwise ENV.sh makes a new one.

set -eu

usage() {
    cat <<EOF
Usage: sbatch runOncoanalyser.sh INPUT.csv [WORK_DIR]

Run nf-core/oncoanalyser 3.0.0 as a SLURM head job.
Submit from the project directory.

  INPUT.csv  oncoanalyser samplesheet (e.g. oncoanalyser_input.csv)
  WORK_DIR   existing /scratch/core001/bic/socci/Oncoanalyser/<YYYYMMDD_HHMMSS>
             folder to -resume; omit to start a new run in a new folder
  -h         show this help
EOF
}

if [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ]; then
    usage
    exit 0
fi

if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    usage >&2
    exit 1
fi

INPUT=$(realpath "$1")
if [ ! -f "$INPUT" ]; then
    echo "INPUT.csv not found: $1" >&2
    exit 1
fi

source ENV.sh
if [ -n "${2:-}" ]; then
    if [ ! -d "$2" ]; then
        echo "WORK_DIR not found: $2" >&2
        exit 1
    fi
    export NXF_WORK=$2
fi

# Never write to the node's /tmp: head job temp files (Nextflow JVM,
# Singularity image pulls) go to node-local /localscratch.
# A new folder per run (mktemp), so this works both via sbatch and
# when run directly (./runOncoanalyser.sh), where SLURM_JOB_ID is
# unset or belongs to the interactive session.
# No automatic cleanup; the folder is left for manual removal.
mkdir -p /localscratch/soccin/oncoanalyser
NXF_LOCAL_TMP=$(mktemp -d -p /localscratch/soccin/oncoanalyser head_XXXXXX)
export TMPDIR=$NXF_LOCAL_TMP
export SINGULARITY_TMPDIR=$NXF_LOCAL_TMP
export APPTAINER_TMPDIR=$NXF_LOCAL_TMP
export NXF_OPTS="$NXF_OPTS -Djava.io.tmpdir=$NXF_LOCAL_TMP"

echo "INPUT=$INPUT"
echo "NXF_WORK=$NXF_WORK"
echo "TMPDIR=$TMPDIR"

./bin/nextflow run nf-core/oncoanalyser \
    -revision 3.0.0 \
    -profile singularity \
    -config config/genome_GRCh38_hmf.config \
    -c config/iris.config \
    -params-file params.yaml \
    --input "$INPUT" \
    -resume
```

`ENV.sh` sets `NXF_WORK`, so no `-work-dir` is needed; `work/` goes to
`/scratch/core001/bic/socci/Oncoanalyser/<YYYYMMDD_HHMMSS>`, a new folder
each time `ENV.sh` is sourced (e.g. `20261112_193000`). `-resume` only
reuses tasks from the same work folder, so to resume pass that folder as
the second argument. The samplesheet and folder in use are printed in
`oncoanalyser.<jobid>.out` as `INPUT=...` and `NXF_WORK=...`; the folder
is also shown as `workDir` in the run header and in `.nextflow.log`.
Results go to `./out`.

The head job writes its temporary files (Nextflow JVM, Singularity image
pulls) to a new `/localscratch/soccin/oncoanalyser/head_XXXXXX` folder (made by
`mktemp`, so it works with `sbatch` or run directly) on its node, never
`/tmp`; the path is printed as `TMPDIR=...` in the job log. No automatic
cleanup.

`-t 7-00:00:00` is the `cmobic_cpu` walltime limit; the head job
must outlive the whole run. If a run needs longer, resubmit with the
work folder as the second argument so `-resume` reuses completed tasks;
resubmitting without it starts a new folder and reruns everything.

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

1. Install Nextflow 25.10.4 in `./bin` (section 2b), set the
   environment variables (section 2c), and check Singularity or
   Apptainer is available on the compute nodes.
2. Inspect the WGS inputs: FASTQ or BAM? If BAM, check reference
   contigs, `-Y`, MC tags, and HLA ALT contigs.
3. Decide GRCh37_hmf vs GRCh38_hmf.
4. Stage reference data (step 4) and write
   `config/genome_GRCh38_hmf.config`.
5. Check `config/iris.config`, write `params.yaml` and the
   samplesheet.
6. Run one tumor/normal pair first; check the ORANGE report and
   PURPLE purity/ploidy before submitting the full cohort.
7. Delete the `work/` directory once outputs are verified.

## 10. Open questions for this project

- Input format of the WGS data (FASTQ vs BAM/CRAM), and the reference
  used if BAM.
- Tumor/normal pairs or tumor-only?
- Is there RNA-seq to include (would switch to WGTS with RNA rows)?
- Target genome build for downstream comparison with other results.
- UMIs in the library prep?
