#!/bin/bash

# Prism: collect the human-readable reports from an oncoanalyser output
# folder into a delivery folder, sorted by file type. Run it from the
# project directory:
#   ./Prism/deliverLocal.sh [OUT_DIR] [RESULTS_DIR]
# OUT_DIR defaults to out (outdir in conf/params.yaml), RESULTS_DIR to
# results. From each group folder in OUT_DIR it takes:
#   PDF/   orange/<sample>.orange.pdf        ORANGE report
#   HTML/  linx/<sample>_linx.html           LINX SV report (self-contained)
#   PNG/   purple/plot/<sample>.circos.png   PURPLE circos plot
#   PNG/   cuppa/<sample>.cuppa.vis.png      CUPPA plot (tumor/normal only)
# File names start with the sample ID, so each type folder is flat.
# Files are hard links, not copies: they take no extra space, but both
# folders must be on the same file system, and a file edited in one
# folder changes in the other. Deleting either name keeps the other.
# Run it again after a rerun to add new files. Links already in place
# are skipped; a different file under an existing name is reported and
# left alone, and the script exits 1.

set -eu

usage() {
    cat <<EOF
Usage: ./Prism/deliverLocal.sh [OUT_DIR] [RESULTS_DIR]

Hard link the ORANGE PDFs, LINX HTML reports, PURPLE circos plots and
CUPPA plots from an oncoanalyser output folder into RESULTS_DIR/PDF,
RESULTS_DIR/HTML and RESULTS_DIR/PNG. Start it from the project
directory.

  OUT_DIR      oncoanalyser output folder (default: out)
  RESULTS_DIR  delivery folder, created if missing (default: results)
  -h           show this help
EOF
}

if [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ]; then
    usage
    exit 0
fi

if [ $# -gt 2 ]; then
    usage >&2
    exit 1
fi

OUT_DIR=${1:-out}
RESULTS_DIR=${2:-results}

if [ ! -d "$OUT_DIR" ]; then
    echo "OUT_DIR not found: $OUT_DIR" >&2
    exit 1
fi

# Patterns that match nothing expand to nothing. Hidden files are not
# matched, so macOS ._* files (written when a share is opened from a
# Mac) are left out.
shopt -s nullglob

STATUS=0

# link_files TYPE LABEL FILE...
# Hard link each FILE into RESULTS_DIR/TYPE and print a count.
link_files() {
    local type=$1 label=$2
    shift 2
    local dest=$RESULTS_DIR/$type
    local src target new=0 old=0
    if [ $# -gt 0 ]; then
        mkdir -p "$dest"
    fi
    for src in "$@"; do
        target=$dest/$(basename "$src")
        if [ -e "$target" ]; then
            if [ "$src" -ef "$target" ]; then
                old=$((old + 1))
            else
                echo "Not replaced, a different file has this name: $target" >&2
                STATUS=1
            fi
            continue
        fi
        # -L: if src is a symlink, link the file it points to.
        if ! ln -L "$src" "$target"; then
            echo "Hard link failed: $src" >&2
            echo "OUT_DIR and RESULTS_DIR must be on the same file system." >&2
            exit 1
        fi
        new=$((new + 1))
    done
    echo "$type/ $label: $# found, $new linked, $old already linked"
}

link_files PDF  "ORANGE report"      "$OUT_DIR"/*/orange/*.orange.pdf
link_files HTML "LINX report"        "$OUT_DIR"/*/linx/*_linx.html
link_files PNG  "PURPLE circos plot" "$OUT_DIR"/*/purple/plot/*.circos.png
link_files PNG  "CUPPA plot"         "$OUT_DIR"/*/cuppa/*.cuppa.vis.png

echo "Results in $RESULTS_DIR"
exit $STATUS
