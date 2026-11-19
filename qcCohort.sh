#!/bin/bash
#SBATCH -J qcCohort
#SBATCH -p cpushort,cmobic_short
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -o qcCohort.%j.out

# Prism: cohort QC with QSEE (Hartwig's QC tool) in multi-sample mode,
# from a finished oncoanalyser run. oncoanalyser runs QSEE once per
# group; this script runs it once over all groups, so the QC status and
# metric values of every sample are in one table and the plot pages in
# one PDF. Each page still shows one tumor (and its normal) against the
# Hartwig cohort percentiles; QSEE 1.0 does not plot samples together.
# Run it from the project directory, as a job or on a compute node:
#   sbatch ./Prism/qcCohort.sh [OUT_DIR]
#   ./Prism/qcCohort.sh [OUT_DIR]
# OUT_DIR defaults to out. ROOT is the folder that holds OUT_DIR; the
# deliverables made after the pipeline go under ROOT/res/prism, this
# script's to ROOT/res/prism/qc/hmf:
#   multisample.qsee.vis.report.<set>.pdf   one page per tumor
#   sample_ids.<set>.tsv                    the sample list given to QSEE
#   data/multisample.qsee.status.<set>.tsv.gz    PASS/WARN/FAIL per metric
#   data/multisample.qsee.vis.data.<set>.tsv.gz  metric values
# Logs go to OUT_DIR/logs/other, with the pipeline's own:
#   qsee_multisample.<set>.log   QSEE log
#   qcCohort.<jobid>.out         this job's SLURM log, moved there at the
#                                end of the job (sbatch only)
# <set> is tumor_normal or tumor_only. QSEE pairs tumors and normals by
# position in the sample list and needs a normal for every tumor or for
# none, so groups with and without a normal are run separately. A set
# with one sample gets QSEE's single-sample file names instead.
# Samples are read from OUT_DIR: the tumor from
# <group>/purple/<tumor>.purple.purity.tsv, the normal from the other
# folder in <group>/bamtools/. The input files are copied (a few MB) to
# a temp folder under /localscratch, which is not removed.

set -eu
set -o pipefail

usage() {
    cat <<EOF
Usage: sbatch ./Prism/qcCohort.sh [OUT_DIR]
       ./Prism/qcCohort.sh [OUT_DIR]

Run QSEE in multi-sample mode over all groups of an oncoanalyser output
folder. Results go to ROOT/res/prism/qc/hmf, where ROOT is the folder
that holds OUT_DIR; logs go to OUT_DIR/logs/other. Start it from the
project directory; run it directly only on a compute node.

  OUT_DIR  oncoanalyser output folder (default: out)
  -h       show this help
EOF
}

if [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ]; then
    usage
    exit 0
fi

if [ $# -gt 1 ]; then
    usage >&2
    exit 1
fi

OUT_DIR=${1:-out}
OUT_DIR=${OUT_DIR%/}
if [ ! -d "$OUT_DIR" ]; then
    echo "OUT_DIR not found: $OUT_DIR" >&2
    exit 1
fi
ROOT=$(cd "$(dirname "$OUT_DIR")" && pwd)
OUT_DIR=$ROOT/$(basename "$OUT_DIR")
QC_DIR=$ROOT/res/prism/qc/hmf
DATA_DIR=$QC_DIR/data
LOG_DIR=$OUT_DIR/logs/other

# Find the Prism folder (PDIR), as in runOncoanalyser.sh: sbatch runs a
# copy of this script from the SLURM spool folder, so $0 only works when
# the script is run directly. The first candidate that holds
# conf/iris.config is used:
#   1. the folder of $0 (direct run)
#   2. SBATCH_SCRIPT_DIR (set by the ~/bin/sbatch wrapper)
#   3. the folder of the submitted script, from scontrol (plain sbatch)
is_prism_dir() {
    [ -n "$1" ] && [ -f "$1/conf/iris.config" ]
}

PDIR=$(dirname "$(readlink -f "$0")")
if ! is_prism_dir "$PDIR"; then
    PDIR=${SBATCH_SCRIPT_DIR:-}
fi
if ! is_prism_dir "$PDIR" && [ -n "${SLURM_JOB_ID:-}" ]; then
    SUBMITTED=$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null \
        | sed -n 's/^ *Command=\([^ ]*\).*/\1/p' || true)
    PDIR=""
    if [ -n "$SUBMITTED" ]; then
        PDIR=$(dirname "$(readlink -f "$SUBMITTED")")
    fi
fi
if ! is_prism_dir "$PDIR"; then
    echo "Cannot find the Prism folder (no conf/iris.config)." >&2
    echo "Set SBATCH_SCRIPT_DIR to the Prism folder and resubmit." >&2
    exit 1
fi

source "$PDIR/conf/ENV.sh"

# Same image as the QSEE module of oncoanalyser 3.0.0
# (modules/local/qsee/main.nf); the pipeline run puts it in the cache.
QSEE_IMG=$NXF_SINGULARITY_CACHEDIR/depot.galaxyproject.org-singularity-hmftools-qsee-1.0--hdfd78af_0.img
if [ ! -f "$QSEE_IMG" ]; then
    echo "QSEE image not found: $QSEE_IMG" >&2
    exit 1
fi

shopt -s nullglob

# Reference data: the HMF resource folder the run used, from the newest
# params file the pipeline wrote. Prism runs GRCh38_hmf only, hence the
# .38 driver gene panel. Relative paths as in oncoanalyser's
# conf/hmf_data.config.
PARAMS_FILES=("$OUT_DIR"/pipeline_info/params_*.json)
if [ ${#PARAMS_FILES[@]} -eq 0 ]; then
    echo "No params_*.json in $OUT_DIR/pipeline_info" >&2
    exit 1
fi
PARAMS=${PARAMS_FILES[${#PARAMS_FILES[@]}-1]}
HMF_DATA=$(sed -n 's/^ *"ref_data_hmf_data_path": *"\([^"]*\)".*$/\1/p' "$PARAMS")
HMF_DATA=${HMF_DATA%/}
DRIVER_GENE_PANEL=$HMF_DATA/common/DriverGenePanel.38.tsv
COHORT_PERCENTILES=$HMF_DATA/misc/qsee/qsee.cohort.percentiles.tsv.gz
for f in "$DRIVER_GENE_PANEL" "$COHORT_PERCENTILES"; do
    if [ ! -f "$f" ]; then
        echo "Reference file not found: $f (from $PARAMS)" >&2
        exit 1
    fi
done

# Never write to the node's /tmp. A new folder per run (mktemp), as in
# runOncoanalyser.sh; it is mounted over the container's /tmp and at its
# own path. No automatic cleanup.
mkdir -p /localscratch/soccin/oncoanalyser
TMP_DIR=$(mktemp -d -p /localscratch/soccin/oncoanalyser qcCohort_XXXXXX)
export TMPDIR=$TMP_DIR
export SINGULARITY_TMPDIR=$TMP_DIR
export APPTAINER_TMPDIR=$TMP_DIR
STAGE_DIR=$TMP_DIR/inputs
mkdir -p "$STAGE_DIR" "$DATA_DIR" "$LOG_DIR"

THREADS=${SLURM_CPUS_PER_TASK:-4}

PRISM_VERSION=$(git -C "$PDIR" describe --always --dirty 2>/dev/null || echo unknown)

echo "PRISM=$PDIR ($PRISM_VERSION)"
echo "OUT_DIR=$OUT_DIR"
echo "QC_DIR=$QC_DIR"
echo "LOG_DIR=$LOG_DIR"
echo "QSEE_IMG=$QSEE_IMG"
echo "HMF_DATA=$HMF_DATA"
echo "TMPDIR=$TMPDIR"

# Find the groups. A group folder has bamtools/; folders without it
# (pipeline_info, logs) are not groups.
TN_GROUPS=()
TN_TUMORS=()
TN_NORMALS=()
TO_GROUPS=()
TO_TUMORS=()
for group_dir in "$OUT_DIR"/*/; do
    group_dir=${group_dir%/}
    if [ ! -d "$group_dir/bamtools" ]; then
        continue
    fi
    purity=("$group_dir"/purple/*.purple.purity.tsv)
    if [ ${#purity[@]} -ne 1 ]; then
        echo "Skipped $group_dir: ${#purity[@]} PURPLE purity files, need 1" >&2
        continue
    fi
    tumor=$(basename "${purity[0]}" .purple.purity.tsv)
    normals=()
    for sample_dir in "$group_dir"/bamtools/*/; do
        sample=$(basename "$sample_dir")
        if [ "$sample" != "$tumor" ]; then
            normals+=("$sample")
        fi
    done
    case ${#normals[@]} in
        0)
            TO_GROUPS+=("$group_dir")
            TO_TUMORS+=("$tumor")
            ;;
        1)
            TN_GROUPS+=("$group_dir")
            TN_TUMORS+=("$tumor")
            TN_NORMALS+=("${normals[0]}")
            ;;
        *)
            echo "More than one normal in $group_dir/bamtools: ${normals[*]}" >&2
            exit 1
            ;;
    esac
done

# stage_sample SAMPLE GROUP_DIR
# Copy the QSEE 1.0 input files of SAMPLE into STAGE_DIR (list from the
# QSEE 1.0 README). Files a sample does not have (PURPLE and ESVEE for
# a normal) are skipped; QSEE runs with -allow_missing_input.
stage_sample() {
    local sample=$1 group_dir=$2 f n=0
    for f in \
        "$group_dir/purple/$sample.purple.purity.tsv" \
        "$group_dir/purple/$sample.purple.qc" \
        "$group_dir/bamtools/$sample/$sample.bam_metric.summary.tsv" \
        "$group_dir/bamtools/$sample/$sample.bam_metric.coverage.tsv" \
        "$group_dir/bamtools/$sample/$sample.bam_metric.flag_counts.tsv" \
        "$group_dir/bamtools/$sample/$sample.bam_metric.frag_length.tsv" \
        "$group_dir/bamtools/$sample/$sample.bam_metric.gene_coverage.tsv" \
        "$group_dir/cobalt/$sample.cobalt.gc.median.tsv" \
        "$group_dir/esvee/$sample.esvee.prep.disc_stats.tsv" \
        "$group_dir/alignments/$sample/$sample.redux.bqr.tsv" \
        "$group_dir/alignments/$sample/$sample.redux.duplicate_freq.tsv" \
        "$group_dir/alignments/$sample/$sample.redux.ms_table.tsv.gz"; do
        if [ -f "$f" ]; then
            cp -L "$f" "$STAGE_DIR/"
            n=$((n + 1))
        fi
    done
    echo "  $sample: $n of 12 input files"
}

# run_qsee SET SAMPLE_ID_FILE
# QSEE writes all its files to one folder and its plot step reads the
# metric table back from there, so it writes to DATA_DIR; the report is
# then moved up to QC_DIR.
run_qsee() {
    local set=$1 ids=$2 f
    echo "Running QSEE: $set"
    singularity exec --no-home --cleanenv \
        --env XDG_CACHE_HOME=/tmp \
        -B "$TMP_DIR:/tmp" -B "$TMP_DIR" -B "$QC_DIR" -B "$HMF_DATA" \
        "$QSEE_IMG" \
        qsee -Xmx12g \
            -sample_id_file "$ids" \
            -sample_data_dir "$STAGE_DIR" \
            -sequencing_type ILLUMINA \
            -driver_gene_panel "$DRIVER_GENE_PANEL" \
            -cohort_percentiles_file "$COHORT_PERCENTILES" \
            -allow_missing_input \
            -merge_plots \
            -threads "$THREADS" \
            -output_id "$set" \
            -output_dir "$DATA_DIR" \
        2>&1 | tee "$LOG_DIR/qsee_multisample.$set.log"
    # The report: multisample.qsee.vis.report.<set>.pdf, or for a set of
    # one sample <tumor>.qsee.vis.report.<set>.pdf and .png.
    for f in "$DATA_DIR"/*.qsee.vis.report."$set".*; do
        mv -f "$f" "$QC_DIR/"
    done
}

if [ ${#TN_TUMORS[@]} -eq 0 ] && [ ${#TO_TUMORS[@]} -eq 0 ]; then
    echo "No groups found in $OUT_DIR" >&2
    exit 1
fi

if [ ${#TN_TUMORS[@]} -gt 0 ]; then
    echo "Tumor/normal groups: ${#TN_TUMORS[@]}"
    ids=$QC_DIR/sample_ids.tumor_normal.tsv
    printf 'TumorId\tReferenceId\n' > "$ids"
    for i in "${!TN_TUMORS[@]}"; do
        printf '%s\t%s\n' "${TN_TUMORS[$i]}" "${TN_NORMALS[$i]}" >> "$ids"
        stage_sample "${TN_TUMORS[$i]}" "${TN_GROUPS[$i]}"
        stage_sample "${TN_NORMALS[$i]}" "${TN_GROUPS[$i]}"
    done
fi

if [ ${#TO_TUMORS[@]} -gt 0 ]; then
    echo "Tumor-only groups: ${#TO_TUMORS[@]}"
    ids=$QC_DIR/sample_ids.tumor_only.tsv
    printf 'TumorId\n' > "$ids"
    for i in "${!TO_TUMORS[@]}"; do
        printf '%s\n' "${TO_TUMORS[$i]}" >> "$ids"
        stage_sample "${TO_TUMORS[$i]}" "${TO_GROUPS[$i]}"
    done
fi

if [ ${#TN_TUMORS[@]} -gt 0 ]; then
    run_qsee tumor_normal "$QC_DIR/sample_ids.tumor_normal.tsv"
fi
if [ ${#TO_TUMORS[@]} -gt 0 ]; then
    run_qsee tumor_only "$QC_DIR/sample_ids.tumor_only.tsv"
fi

echo "Results in $QC_DIR"
echo "Logs in $LOG_DIR"

# Move this job's SLURM log (#SBATCH -o, in the submit folder) to
# LOG_DIR. Only a file with exactly that name is moved, so a direct run
# inside an interactive allocation leaves that session's log alone. On
# the same file system mv is a rename and SLURM keeps writing to the
# moved file, so nothing after this line is lost.
SLURM_LOG=qcCohort.${SLURM_JOB_ID:-none}.out
if [ -n "${SLURM_JOB_ID:-}" ] && [ -f "$SLURM_LOG" ]; then
    echo "SLURM log moved to $LOG_DIR/$SLURM_LOG"
    mv "$SLURM_LOG" "$LOG_DIR/"
fi
