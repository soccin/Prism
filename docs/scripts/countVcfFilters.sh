#!/bin/bash
#SBATCH -J countVcfFilters
#SBATCH -p cpushort,cmobic_short
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 02:00:00
#SBATCH -o countVcfFilters.%j.out

# Prism: count the records of every variant VCF along each calling chain,
# by FILTER value, for every group of a finished oncoanalyser run. Used
# to confirm what each step keeps and drops (docs/ACTUAL_METHODS.md).
# Read-only on OUT_DIR. Run from the project directory:
#   sbatch ./Prism/docs/scripts/countVcfFilters.sh OUT_DIR RESULT_DIR
# Writes to RESULT_DIR:
#   filter_counts.tsv   group, mode, step, file, filter, records
#   svtype_counts.tsv   group, mode, step, file, svtype, filter, records
#                       (SV files only)
# mode is matched (a germline PURPLE VCF exists) or tumor_only. step names
# the chain position, e.g. snv_somatic.1_sage. A FILTER value with
# several labels (a;b) is counted as written. Missing files are skipped.
# bcftools 1.20 is taken from the module path if not on PATH. No temp
# files are written.

set -eu
set -o pipefail

if [ $# -ne 2 ]; then
    echo "Usage: sbatch ./Prism/docs/scripts/countVcfFilters.sh OUT_DIR RESULT_DIR" >&2
    exit 1
fi
OUT_DIR=$1
RESULT_DIR=$2

if ! command -v bcftools > /dev/null; then
    export PATH=/admin/software/bcftools/bcftools-1.20/bin:$PATH
fi
echo "bcftools: $(command -v bcftools) $(bcftools --version | head -1)"

mkdir -p "$RESULT_DIR"
FILTERS=$RESULT_DIR/filter_counts.tsv
SVTYPES=$RESULT_DIR/svtype_counts.tsv
printf "group\tmode\tstep\tfile\tfilter\trecords\n" > "$FILTERS"
printf "group\tmode\tstep\tfile\tsvtype\tfilter\trecords\n" > "$SVTYPES"

# step and path (relative to the group folder, <T> = tumor ID)
CHAINS="
snv_somatic.1_sage      sage/somatic/<T>.sage.somatic.vcf.gz
snv_somatic.2_pave      pave/somatic/<T>.pave.somatic.vcf.gz
snv_somatic.3_purple    purple/<T>.purple.somatic.vcf.gz
snv_germline.1_sage     sage/germline/<T>.sage.germline.vcf.gz
snv_germline.2_pave     pave/germline/<T>.pave.germline.vcf.gz
snv_germline.3_purple   purple/<T>.purple.germline.vcf.gz
sv.1_esvee_raw          esvee/<T>.esvee.raw.vcf.gz
sv.2_esvee_unfiltered   esvee/<T>.esvee.unfiltered.vcf.gz
sv_somatic.3_esvee      esvee/<T>.esvee.somatic.vcf.gz
sv_somatic.4_purple     purple/<T>.purple.sv.vcf.gz
sv_germline.3_esvee     esvee/<T>.esvee.germline.vcf.gz
sv_germline.4_purple    purple/<T>.purple.sv.germline.vcf.gz
"

for purity in "$OUT_DIR"/*/purple/*.purple.purity.tsv; do
    group_dir=$(dirname "$(dirname "$purity")")
    group=$(basename "$group_dir")
    tumor=$(basename "$purity" .purple.purity.tsv)
    if [ -e "$group_dir/purple/$tumor.purple.germline.vcf.gz" ]; then
        mode=matched
    else
        mode=tumor_only
    fi
    echo "$(date '+%H:%M:%S') $group ($mode)"
    while read -r step pattern; do
        [ -z "$step" ] && continue
        rel=${pattern//<T>/$tumor}
        file=$group_dir/$rel
        [ -e "$file" ] || continue
        name=${rel//$tumor/<T>}
        bcftools query -f '%FILTER\n' "$file" \
            | sort | uniq -c \
            | awk -v g="$group" -v m="$mode" -v s="$step" -v f="$name" \
                'BEGIN { OFS = "\t" } { print g, m, s, f, $2, $1 }' >> "$FILTERS"
        if [[ $step == sv* ]]; then
            bcftools query -f '%INFO/SVTYPE\t%FILTER\n' "$file" \
                | sort | uniq -c \
                | awk -v g="$group" -v m="$mode" -v s="$step" -v f="$name" \
                    'BEGIN { OFS = "\t" } { print g, m, s, f, $2, $3, $1 }' >> "$SVTYPES"
        fi
    done <<< "$CHAINS"
done

echo "$(date '+%H:%M:%S') done: $FILTERS $SVTYPES"
