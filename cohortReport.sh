#!/bin/bash
#SBATCH -J cohortReport
#SBATCH -p cpushort,cmobic_short
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -o cohortReport.%j.out

# Prism: cohort report from a finished oncoanalyser run. Reads the
# per-sample result tables of every group and writes one workbook with
# the events per sample and their counts over the cohort: somatic
# variants, gene copy number, gene disruptions with breakpoints and
# support, fusions, viruses, whole genome duplication, DPYD status.
# It counts and lists; it applies no significance test and leaves no
# sample out. Run it from the project directory, as a job or on a
# compute node:
#   sbatch ./Prism/cohortReport.sh [-m MANIFEST] [-g GENES] [-p FRAC] [OUT_DIR]
#   ./Prism/cohortReport.sh [-m MANIFEST] [-g GENES] [-p FRAC] [OUT_DIR]
# OUT_DIR defaults to out. ROOT is the folder that holds OUT_DIR; the
# report goes to ROOT/res/prism/cohort:
#   prism_report.xlsx   one sheet per table, the genes of interest
#                       first; sheet "manifest" says what each sheet
#                       holds, sheet "glossary" explains every column
#   prism_sv_same_breakpoint.xlsx
#                       breakpoints found at the same base pair in more
#                       than one sample, for review
#   genes_of_interest_default.txt
#                       the default genes of interest (only without -g)
#   data/<table>.csv    the same tables as CSV, and gene_summary_all.csv
#                       with every gene
# The gene tables in the workbook keep a gene altered in at least
# max(2, FRAC x samples) samples (-p, default 0.05). gene_summary_HC
# is counted without the breakends at a breakpoint that more than one
# sample has; gene_summary_UnFilt counts every breakend.
# Logs go to OUT_DIR/logs/other, with the pipeline's own:
#   cohortReport.<jobid>.out   this job's SLURM log, moved there at the
#                              end of the job (sbatch only)
# MANIFEST is a CSV with the columns sample_id (tumor sample) and
# sample_class (for example Tumor, CellLine); other columns are copied
# to the samples sheet. Without it every sample has class Tumor. Events
# on chrX and chrY are reported only for samples with a value in a
# manifest column named sex and a matched normal.
# GENES is a text file with one gene symbol per line (# starts a
# comment): the genes of interest. They are marked in the gene tables and
# get the first two sheets: genes_of_interest_events (every event, one
# row per event) and genes_of_interest (every gene in every sample).
# Without -g the genes of interest are the genes with a pipeline driver
# call of likelihood 0.8 or higher in at least one sample; that list is
# written to genes_of_interest_default.txt, to be edited and given back
# with -g.
# A rerun replaces the two workbooks, the CSV files and, without -g, the
# default gene list.

set -eu
set -o pipefail

usage() {
    cat <<EOF
Usage: sbatch ./Prism/cohortReport.sh [-m MANIFEST] [-g GENES] [-p FRAC] [OUT_DIR]
       ./Prism/cohortReport.sh [-m MANIFEST] [-g GENES] [-p FRAC] [OUT_DIR]

Write the cohort report of an oncoanalyser output folder to
ROOT/res/prism/cohort (prism_report.xlsx, prism_sv_same_breakpoint.xlsx
and data/*.csv), where ROOT is the folder that holds OUT_DIR. Start it
from the project directory; run it directly only on a compute node.

  -m MANIFEST  CSV with sample_id and sample_class; optional column sex
  -g GENES     text file, one gene of interest per line (default: the
               genes with a pipeline driver call of likelihood 0.8 or
               higher; written to genes_of_interest_default.txt)
  -p FRAC      gene tables keep a gene altered in at least
               max(2, FRAC x samples) samples (default: 0.05)
  OUT_DIR      oncoanalyser output folder (default: out)
  -h           show this help
EOF
}

MANIFEST=""
GENES=""
MIN_PCT=""
while getopts ":m:g:p:h" opt; do
    case $opt in
        m) MANIFEST=$OPTARG ;;
        g) GENES=$OPTARG ;;
        p) MIN_PCT=$OPTARG ;;
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
RES_DIR=$ROOT/res/prism/cohort
LOG_DIR=$OUT_DIR/logs/other

R_ARGS=()
if [ -n "$MANIFEST" ]; then
    if [ ! -f "$MANIFEST" ]; then
        echo "MANIFEST not found: $MANIFEST" >&2
        exit 1
    fi
    R_ARGS+=(--manifest "$(readlink -f "$MANIFEST")")
fi
if [ -n "$GENES" ]; then
    if [ ! -f "$GENES" ]; then
        echo "GENES not found: $GENES" >&2
        exit 1
    fi
    R_ARGS+=(--genes "$(readlink -f "$GENES")")
fi
if [ -n "$MIN_PCT" ]; then
    R_ARGS+=(--min-pct "$MIN_PCT")
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

shopt -s nullglob

# The driver gene panel of the HMF resource folder the run used, from
# the newest params file the pipeline wrote. Prism runs GRCh38_hmf only,
# hence the .38 panel. Without it the report has no panel columns.
PARAMS_FILES=("$OUT_DIR"/pipeline_info/params_*.json)
if [ ${#PARAMS_FILES[@]} -gt 0 ]; then
    PARAMS=${PARAMS_FILES[${#PARAMS_FILES[@]}-1]}
    HMF_DATA=$(sed -n 's/^ *"ref_data_hmf_data_path": *"\([^"]*\)".*$/\1/p' "$PARAMS")
    PANEL=${HMF_DATA%/}/common/DriverGenePanel.38.tsv
    if [ -f "$PANEL" ]; then
        R_ARGS+=(--panel "$PANEL")
    else
        echo "Driver gene panel not found: $PANEL (from $PARAMS)" >&2
    fi
else
    echo "No params_*.json in $OUT_DIR/pipeline_info: no driver gene panel" >&2
fi

# Never write to the node's /tmp. A new folder per run (mktemp), as in
# runOncoanalyser.sh. No automatic cleanup.
mkdir -p /localscratch/soccin/oncoanalyser
TMP_DIR=$(mktemp -d -p /localscratch/soccin/oncoanalyser cohortReport_XXXXXX)
export TMPDIR=$TMP_DIR

mkdir -p "$RES_DIR/data" "$LOG_DIR"

PRISM_VERSION=$(git -C "$PDIR" describe --always --dirty 2>/dev/null || echo unknown)

echo "PRISM=$PDIR ($PRISM_VERSION)"
echo "OUT_DIR=$OUT_DIR"
echo "RES_DIR=$RES_DIR"
echo "LOG_DIR=$LOG_DIR"
echo "MANIFEST=${MANIFEST:-none}"
echo "GENES=${GENES:-none (default list)}"
echo "TMPDIR=$TMPDIR"

Rscript --vanilla "$PDIR/scripts/cohortReport.R" \
    --out-dir "$OUT_DIR" \
    --res-dir "$RES_DIR" \
    --version "$PRISM_VERSION" \
    ${R_ARGS[@]+"${R_ARGS[@]}"}

echo "Report in $RES_DIR"
echo "Logs in $LOG_DIR"

# Move this job's SLURM log (#SBATCH -o, in the submit folder) to
# LOG_DIR. Only a file with exactly that name is moved, so a direct run
# inside an interactive allocation leaves that session's log alone. On
# the same file system mv is a rename and SLURM keeps writing to the
# moved file, so nothing after this line is lost.
SLURM_LOG=cohortReport.${SLURM_JOB_ID:-none}.out
if [ -n "${SLURM_JOB_ID:-}" ] && [ -f "$SLURM_LOG" ]; then
    echo "SLURM log moved to $LOG_DIR/$SLURM_LOG"
    mv "$SLURM_LOG" "$LOG_DIR/"
fi
