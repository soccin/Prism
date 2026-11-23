#!/bin/bash
#SBATCH -J postOA
#SBATCH -p cpushort,cmobic_short
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH -o postOA.%j.out

# Prism: cohort tables of a finished oncoanalyser run, for analysis over
# the cohort. Four steps:
#   1. scripts/postOA.R stacks the per-sample LINX tables of every group
#      into one table each, with the group and the sample in front. Lines
#      are copied as LINX wrote them: no row is left out and no value is
#      changed.
#   2. scripts/svCohort.R writes one table of the somatic SV events of all
#      groups: every record of PURPLE's purple.sv.vcf.gz (all FILTER
#      values), one row per event, with every field PURPLE, ESVEE and LINX
#      wrote for it, region tracks (TRACK_DIR) and the same breakends in
#      other groups, unpacked for analysis (fixed pairs split, repeated
#      columns dropped; each step checked). No row is left out and no
#      value is interpreted.
#   3. scripts/tsvToParquet.py (python3 with pyarrow) writes the typed
#      Parquet copy of that table, with the types in its column dictionary.
#   4. docs/scripts/checkSvCohort.R rebuilds every PURPLE VCF record and
#      ESVEE row from the table and checks the Parquet copy; the job fails
#      if any check fails.
# Run it from the project directory, as a job or on a compute node:
#   sbatch ./Prism/postOA.sh [-t TRACK_DIR] [OUT_DIR]
#   ./Prism/postOA.sh [-t TRACK_DIR] [OUT_DIR]
# OUT_DIR defaults to out. TRACK_DIR is the folder written by
# scripts/fetchTracks.sh (default below). ROOT is the folder that holds
# OUT_DIR; the tables go to ROOT/res/prism/db:
#   manifest.tsv         what each file holds, where it comes from, and
#                        the number of groups and rows in it
#   samples.tsv          one row per group: group_id, sample_id (tumor),
#                        normal_id, run_mode, linx_version
#   rows_per_sample.tsv  the number of rows each group adds to each table
#   sv_cohort_filtered.parquet
#                        one row per somatic SV event (steps 2, 3), typed;
#                        lists as arrays. Read with nanoparquet (R),
#                        pyarrow or polars (Python)
#   sv_cohort_filtered.tsv.gz
#                        the same table as text; lists joined with | (VCF
#                        lists with ,)
#   sv_cohort_filtered.columns.tsv
#                        its column dictionary: description, type, list
#                        separator and count column, and the number of
#                        non-empty values of every column
#   data/linx_somatic_<table>.tsv.gz
#                        every group: svs, clusters, links, breakend,
#                        fusion, drivers, driver_catalog
#   data/linx_germline_<table>.tsv.gz
#                        groups with a matched normal: svs, breakend,
#                        disruption, clusters, links, driver_catalog
# Each data table has the columns group_id and sample_id, then the
# columns of the LINX file (see the LINX README). LINX's neoepitope and
# vis_* files are not aggregated. The job stops if a file is missing, if
# the files of a table have different headers, if a line has a different
# number of fields than its header, if the SV events do not account for
# every PURPLE record and LINX row, if an unpacking step finds a row that
# breaks its assumption, if a value does not fit its column's type, or if
# a check of step 4 fails.
# The Hartwig resource folder (sv_prep_blacklist, unmap_regions) is read
# from ref_data_hmf_data_path in the newest OUT_DIR/pipeline_info/params_*.json.
# Logs go to OUT_DIR/logs/other, with the pipeline's own:
#   postOA.<jobid>.out   this job's SLURM log, moved there at the end of
#                        the job (sbatch only)
# A rerun replaces these files.

set -eu
set -o pipefail

usage() {
    cat <<EOF
Usage: sbatch ./Prism/postOA.sh [-t TRACK_DIR] [OUT_DIR]
       ./Prism/postOA.sh [-t TRACK_DIR] [OUT_DIR]

Write the cohort tables of an oncoanalyser output folder into
ROOT/res/prism/db, where ROOT is the folder that holds OUT_DIR: the
stacked LINX tables (data/*.tsv.gz) and one table of the somatic SV
events (sv_cohort_filtered.parquet and .tsv.gz), checked against its
sources. Needs Rscript (tidyverse, fs, GenomicRanges, nanoparquet) and
python3 with pyarrow. Start it from the project directory; run it
directly only on a compute node.

  -t TRACK_DIR region tracks from scripts/fetchTracks.sh
               (default: $DEFAULT_TRACK_DIR)
  OUT_DIR      oncoanalyser output folder (default: out)
  -h           show this help
EOF
}

DEFAULT_TRACK_DIR=/scratch/core001/bic/socci/reference/prism_tracks/GRCh38
TRACK_DIR=$DEFAULT_TRACK_DIR
while getopts ":t:h" opt; do
    case $opt in
        t) TRACK_DIR=$OPTARG ;;
        h) usage; exit 0 ;;
        *) usage >&2; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

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
RES_DIR=$ROOT/res/prism/db
LOG_DIR=$OUT_DIR/logs/other

TRACK_DIR=${TRACK_DIR%/}
if [ ! -f "$TRACK_DIR/tracks.tsv" ]; then
    echo "No tracks.tsv in TRACK_DIR: $TRACK_DIR (run scripts/fetchTracks.sh)" >&2
    exit 1
fi

# The Hartwig resource folder of the run, from the newest params file.
shopt -s nullglob
PARAMS_FILES=("$OUT_DIR"/pipeline_info/params_*.json)
shopt -u nullglob
if [ ${#PARAMS_FILES[@]} -eq 0 ]; then
    echo "No params_*.json in $OUT_DIR/pipeline_info" >&2
    exit 1
fi
PARAMS=${PARAMS_FILES[${#PARAMS_FILES[@]}-1]}
HMF_DIR=$(sed -n 's/^ *"ref_data_hmf_data_path": *"\([^"]*\)".*$/\1/p' "$PARAMS")
HMF_DIR=${HMF_DIR%/}
if [ ! -d "$HMF_DIR" ]; then
    echo "Hartwig resource folder not found: '$HMF_DIR' (from $PARAMS)" >&2
    exit 1
fi

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

if ! command -v Rscript >/dev/null; then
    echo "Rscript not found on PATH" >&2
    exit 1
fi
if ! Rscript --vanilla -e 'quit(status = !requireNamespace("nanoparquet", quietly = TRUE))'; then
    echo "R package nanoparquet not installed (needed by checkSvCohort.R)" >&2
    exit 1
fi
if ! python3 -c 'import pyarrow' 2>/dev/null; then
    echo "python3 with pyarrow not found (python3 -m pip install --user pyarrow)" >&2
    exit 1
fi

# Never write to the node's /tmp. A new folder per run (mktemp), as in
# runOncoanalyser.sh. No automatic cleanup.
mkdir -p /localscratch/soccin/oncoanalyser
TMP_DIR=$(mktemp -d -p /localscratch/soccin/oncoanalyser postOA_XXXXXX)
export TMPDIR=$TMP_DIR

mkdir -p "$RES_DIR/data" "$LOG_DIR"

PRISM_VERSION=$(git -C "$PDIR" describe --always --dirty 2>/dev/null || echo unknown)

echo "PRISM=$PDIR ($PRISM_VERSION)"
echo "OUT_DIR=$OUT_DIR"
echo "RES_DIR=$RES_DIR"
echo "LOG_DIR=$LOG_DIR"
echo "TMPDIR=$TMPDIR"
echo "TRACK_DIR=$TRACK_DIR"
echo "HMF_DIR=$HMF_DIR"

Rscript --vanilla "$PDIR/scripts/postOA.R" \
    --out-dir "$OUT_DIR" \
    --res-dir "$RES_DIR"

Rscript --vanilla "$PDIR/scripts/svCohort.R" \
    --out-dir "$OUT_DIR" \
    --db-dir "$RES_DIR" \
    --track-dir "$TRACK_DIR" \
    --hmf-dir "$HMF_DIR"

python3 "$PDIR/scripts/tsvToParquet.py" \
    "$RES_DIR/sv_cohort_filtered.tsv.gz" \
    "$RES_DIR/sv_cohort_filtered.columns.tsv" \
    "$RES_DIR/sv_cohort_filtered.parquet"

Rscript --vanilla "$PDIR/docs/scripts/checkSvCohort.R" "$OUT_DIR" "$RES_DIR"

echo "Tables in $RES_DIR"
echo "Logs in $LOG_DIR"

# Move this job's SLURM log (#SBATCH -o, in the submit folder) to
# LOG_DIR. Only a file with exactly that name is moved, so a direct run
# inside an interactive allocation leaves that session's log alone. On
# the same file system mv is a rename and SLURM keeps writing to the
# moved file, so nothing after this line is lost.
SLURM_LOG=postOA.${SLURM_JOB_ID:-none}.out
if [ -n "${SLURM_JOB_ID:-}" ] && [ -f "$SLURM_LOG" ]; then
    echo "SLURM log moved to $LOG_DIR/$SLURM_LOG"
    mv "$SLURM_LOG" "$LOG_DIR/"
fi
