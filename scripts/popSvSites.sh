#!/bin/bash
#SBATCH -J popSvSites
#SBATCH -p cpushort,cmobic_short
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 01:30:00
#SBATCH -o popSvSites.%j.out

# Prism: reduce the population SV VCFs (short-read WGS call sets,
# downloaded by fetchTracks.sh) to the fields svCohort.R matches SV
# events against. One table per source, one line per VCF record, every
# record kept (any FILTER):
#   id chrom pos end chr2 end2 svtype svlen filter af
# Missing values are "." (as bcftools query writes them).
#   sbatch ./Prism/scripts/popSvSites.sh [TRACK_DIR]
# TRACK_DIR defaults to the one of fetchTracks.sh. Writes
# TRACK_DIR/<source>/sites.tsv.gz for each source below. BCFTOOLS
# (default below) is the bcftools binary.
#
# Sources and their fields:
#   gnomad_sv      gnomAD-SV v4.1: END, CHR2, END2, SVTYPE, SVLEN, FILTER, AF
#   kgp_sv         1000 Genomes high coverage: the same fields
#   dbvar_nstd166  gnomAD-SV v2.1 (dbVar): END, CHR2, SVTYPE, SVLEN, AF;
#   dbvar_nstd200  CCDG (dbVar): the same. dbVar has no END2 and no
#                  FILTER ("."). For a record whose ALT is in breakend
#                  notation (BND, one record per junction), end2 is the
#                  partner position in ALT, and chr2 its chromosome where
#                  CHR2 is missing. dbVar names the primary
#                  chromosomes 1-22, X, Y, MT: these are written chr1-22,
#                  chrX, chrY, chrM, as in GRCh38 of the run. Other
#                  contigs (NT_, NW_, GL accessions) are kept as written;
#                  they match no breakend.
# The job stops if a table's line count differs from the record count of
# the VCF: from its index, or, where the index has no count metadata
# (written by an older tabix, as the 1000 Genomes one), its data lines.

set -eu
set -o pipefail

TRACK_DIR=${1:-/scratch/core001/bic/socci/reference/prism_tracks/GRCh38}
TRACK_DIR=${TRACK_DIR%/}
BCFTOOLS=${BCFTOOLS:-/admin/software/bcftools/bcftools-1.20/bin/bcftools}

if [ ! -x "$BCFTOOLS" ]; then
    echo "bcftools not found: $BCFTOOLS (set BCFTOOLS)" >&2
    exit 1
fi

# Never write to the node's /tmp.
mkdir -p "/localscratch/soccin/popSvSites/${SLURM_JOB_ID:-local}"
export TMPDIR=/localscratch/soccin/popSvSites/${SLURM_JOB_ID:-local}

# Query formats, used by name through ${!FORMAT} below.
# shellcheck disable=SC2034
FULL='%ID\t%CHROM\t%POS\t%INFO/END\t%INFO/CHR2\t%INFO/END2\t%INFO/SVTYPE\t%INFO/SVLEN\t%FILTER\t%INFO/AF\n'
# dbVar: ALT in the end2 column, replaced by the awk below.
# shellcheck disable=SC2034
DBVAR='%ID\t%CHROM\t%POS\t%INFO/END\t%INFO/CHR2\t%ALT\t%INFO/SVTYPE\t%INFO/SVLEN\t%FILTER\t%INFO/AF\n'

# source<TAB>VCF (relative to TRACK_DIR)<TAB>format<TAB>rename chromosomes
SOURCES=$(cat <<EOF
gnomad_sv	gnomad_sv/gnomad.v4.1.sv.sites.vcf.gz	FULL	no
kgp_sv	kgp_sv/1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz	FULL	no
dbvar_nstd166	dbvar_nstd166/nstd166.GRCh38.variant_call.vcf.gz	DBVAR	yes
dbvar_nstd200	dbvar_nstd200/nstd200.GRCh38.variant_call.vcf.gz	DBVAR	yes
EOF
)

"$BCFTOOLS" --version | head -n 1

while IFS=$'\t' read -r NAME FILE FORMAT RENAME; do
    VCF=$TRACK_DIR/$FILE
    OUT=$TRACK_DIR/$NAME/sites.tsv.gz
    if [ ! -s "$VCF" ] || [ ! -s "$VCF.tbi" ]; then
        echo "Not found: $VCF or its .tbi (run scripts/fetchTracks.sh)" >&2
        exit 1
    fi
    echo "$NAME: $VCF"
    FMT=${!FORMAT}
    {
        printf 'id\tchrom\tpos\tend\tchr2\tend2\tsvtype\tsvlen\tfilter\taf\n'
        if [ "$RENAME" = yes ]; then
            "$BCFTOOLS" query -f "$FMT" "$VCF" | awk 'BEGIN {FS = OFS = "\t"}
                function chr(c) { if (c == "MT") return "chrM"; if (c ~ /^([0-9]+|X|Y)$/) return "chr" c; return c }
                {
                    alt = $6; $6 = "."
                    if (match(alt, /[\[\]][^\[\]:]+:[0-9]+[\[\]]/)) {
                        mate = substr(alt, RSTART + 1, RLENGTH - 2)
                        k = split(mate, m, ":")
                        if ($5 == ".") $5 = m[1]
                        $6 = m[2]
                    }
                    $2 = chr($2); if ($5 != ".") $5 = chr($5); print
                }'
        else
            "$BCFTOOLS" query -f "$FMT" "$VCF"
        fi
    } | gzip > "$OUT.part"
    # Record count from the index. An index written by an older tabix has
    # no count metadata (bcftools index -n then prints 0): count the VCF's
    # data lines instead.
    if "$BCFTOOLS" index -s "$VCF" > /dev/null 2>&1; then
        N_VCF=$("$BCFTOOLS" index -n "$VCF")
        HOW="from the VCF index"
    else
        N_VCF=$(zcat "$VCF" | grep -vc '^#')
        HOW="data lines of the VCF; its index has no counts"
    fi
    N_OUT=$(( $(zcat "$OUT.part" | wc -l) - 1 ))
    if [ "$N_VCF" -ne "$N_OUT" ]; then
        echo "$NAME: $N_VCF records ($HOW), but the table has $N_OUT lines" >&2
        exit 1
    fi
    mv "$OUT.part" "$OUT"
    echo "  $N_OUT sites written to $OUT; $N_VCF records ($HOW)"
done <<< "$SOURCES"
