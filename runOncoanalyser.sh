#!/bin/bash
#SBATCH -J oncoanalyser
#SBATCH -p cmobic_cpu
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 7-00:00:00
#SBATCH -o oncoanalyser.%j.out

# Prism: Nextflow head job for nf-core/oncoanalyser 3.0.0.
# Clone Prism into the project directory, run 00.SETUP.sh once, then
# from the project directory either submit it or run it directly:
#   sbatch ./Prism/runOncoanalyser.sh INPUT.csv [WORK_DIR]
#   ./Prism/runOncoanalyser.sh INPUT.csv [WORK_DIR]
# INPUT.csv: the oncoanalyser samplesheet (e.g. oncoanalyser_input.csv).
# Optional WORK_DIR: an existing
#   /scratch/core001/bic/socci/Oncoanalyser/<YYYYMMDD_HHMMSS>
# folder to -resume a previous run; otherwise conf/ENV.sh makes a new one.
# Outputs (./out, oncoanalyser.<jobid>.out, .nextflow.log) go to the
# directory the script is started from, not to the Prism folder.

set -eu

usage() {
    cat <<EOF
Usage: sbatch ./Prism/runOncoanalyser.sh INPUT.csv [WORK_DIR]
       ./Prism/runOncoanalyser.sh INPUT.csv [WORK_DIR]

Run nf-core/oncoanalyser 3.0.0, as a SLURM head job (sbatch) or in
the current shell. Start it from the project directory.

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

# Find the Prism folder (PDIR). sbatch runs a copy of this script from
# the SLURM spool folder, so $0 only works when the script is run
# directly. Candidates, in order; the first one that holds
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

if [ ! -x "$PDIR/bin/nextflow" ]; then
    echo "Nextflow is not installed in $PDIR/bin." >&2
    echo "Run: bash $PDIR/00.SETUP.sh" >&2
    exit 1
fi

source "$PDIR/conf/ENV.sh"
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

PRISM_VERSION=$(git -C "$PDIR" describe --always --dirty 2>/dev/null || echo unknown)

echo "PRISM=$PDIR ($PRISM_VERSION)"
echo "INPUT=$INPUT"
echo "NXF_WORK=$NXF_WORK"
echo "TMPDIR=$TMPDIR"

"$PDIR/bin/nextflow" run nf-core/oncoanalyser \
    -revision 3.0.0 \
    -profile singularity \
    -config "$PDIR/conf/genome_GRCh38_hmf.config" \
    -c "$PDIR/conf/iris.config" \
    -params-file "$PDIR/conf/params.yaml" \
    --input "$INPUT" \
    -resume
