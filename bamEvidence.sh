#!/bin/bash
#SBATCH -J bamEvidence
#SBATCH -p cpushort,cmobic_short
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH -o bamEvidence.%A_%a.out

# Prism: recount the reads at every SV breakend of the cohort in every BAM
# of the run (each tumor and each normal), with scripts/bamEvidence.py.
# One array task per BAM; svCohort.R (postOA.sh) then reads the counts
# into the SV table: the group's own tumor, its own normal, and the other
# samples at the same breakend.
#
# Run it from the project directory after postOA.R has written
# ROOT/res/prism/db/samples.tsv (a first postOA.sh run, or postOA.R
# alone). The array size must equal the number of BAMs; the script prints
# the command when run without SLURM:
#   ./Prism/bamEvidence.sh [OUT_DIR]          print the sbatch command
#   sbatch --array=0-N ./Prism/bamEvidence.sh [OUT_DIR]
# OUT_DIR defaults to out. BAMs: OUT_DIR/<group>/alignments/<sample>/
# <sample>.redux.bam, in sorted order. Output:
#   ROOT/res/prism/db/bam_evidence/<sample>.bam_evidence.tsv.gz
# Logs go to OUT_DIR/logs/other (moved at the end of each task).
# A task stops if its BAM's read groups do not all name the sample, or if
# the array size differs from the number of BAMs.

set -eu
set -o pipefail

OUT_DIR=${1:-out}
OUT_DIR=${OUT_DIR%/}
if [ ! -d "$OUT_DIR" ]; then
    echo "OUT_DIR not found: $OUT_DIR" >&2
    exit 1
fi
ROOT=$(cd "$(dirname "$OUT_DIR")" && pwd)
OUT_DIR=$ROOT/$(basename "$OUT_DIR")
RES_DIR=$ROOT/res/prism/db
EVIDENCE_DIR=$RES_DIR/bam_evidence
LOG_DIR=$OUT_DIR/logs/other
SAMPLES=$RES_DIR/samples.tsv

if [ ! -s "$SAMPLES" ]; then
    echo "No $SAMPLES (run postOA.sh or scripts/postOA.R first)" >&2
    exit 1
fi

shopt -s nullglob
BAMS=("$OUT_DIR"/*/alignments/*/*.redux.bam)
shopt -u nullglob
N=${#BAMS[@]}
if [ "$N" -eq 0 ]; then
    echo "No BAMs under $OUT_DIR/*/alignments" >&2
    exit 1
fi

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "$N BAMs under $OUT_DIR. Submit:"
    echo "  sbatch --array=0-$((N - 1)) $0 $OUT_DIR"
    exit 0
fi
if [ "${SLURM_ARRAY_TASK_COUNT:-0}" -ne "$N" ]; then
    echo "Array has ${SLURM_ARRAY_TASK_COUNT:-0} tasks, but there are $N BAMs" >&2
    exit 1
fi

BAM=${BAMS[$SLURM_ARRAY_TASK_ID]}
SAMPLE=$(basename "$BAM" .redux.bam)

# Find the Prism folder (PDIR), as in postOA.sh.
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

if ! python3 -c 'import pysam' 2>/dev/null; then
    echo "python3 with pysam not found (python3 -m pip install --user pysam)" >&2
    exit 1
fi

# Never write to the node's /tmp.
mkdir -p /localscratch/soccin/oncoanalyser
TMP_DIR=$(mktemp -d -p /localscratch/soccin/oncoanalyser bamEvidence_XXXXXX)
export TMPDIR=$TMP_DIR

mkdir -p "$EVIDENCE_DIR" "$LOG_DIR"

PRISM_VERSION=$(git -C "$PDIR" describe --always --dirty 2>/dev/null || echo unknown)
echo "PRISM=$PDIR ($PRISM_VERSION)"
echo "TASK=$SLURM_ARRAY_TASK_ID of $N"
echo "BAM=$BAM"
echo "SAMPLE=$SAMPLE"
echo "TMPDIR=$TMPDIR"

python3 "$PDIR/scripts/bamEvidence.py" \
    --out-dir "$OUT_DIR" \
    --samples "$SAMPLES" \
    --bam "$BAM" \
    --sample-id "$SAMPLE" \
    --output "$EVIDENCE_DIR/$SAMPLE.bam_evidence.tsv.gz" \
    --threads 2

# Move this task's SLURM log to LOG_DIR (see postOA.sh).
SLURM_LOG=bamEvidence.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out
if [ -f "$SLURM_LOG" ]; then
    echo "SLURM log moved to $LOG_DIR/$SLURM_LOG"
    mv "$SLURM_LOG" "$LOG_DIR/"
fi
