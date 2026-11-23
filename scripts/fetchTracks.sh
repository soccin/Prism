#!/bin/bash
# Prism: download the GRCh38 region tracks that postOA.sh uses to
# annotate SV breakends with technical context (mappability, repeats,
# segmental duplications, blacklists, gaps). Run once per site, on the
# login node (downloads only):
#   ./Prism/scripts/fetchTracks.sh [TRACK_DIR]
# TRACK_DIR defaults to /scratch/core001/bic/socci/reference/prism_tracks/GRCh38.
# Each track goes into its own folder TRACK_DIR/<name>/. TRACK_DIR/tracks.tsv
# lists every track: name, file, source URL, download date, md5, and the
# result of the md5 check against the GIAB md5 list (GIAB files only).
# Files already present are kept and not downloaded again; delete a
# track's folder to fetch it again.

set -eu
set -o pipefail

TRACK_DIR=${1:-/scratch/core001/bic/socci/reference/prism_tracks/GRCh38}
TRACK_DIR=${TRACK_DIR%/}

GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all
UCSC=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database

# name<TAB>URL. The names are the column stems svCohort.R uses.
TRACKS=$(cat <<EOF
giab_alldifficult	$GIAB/Union/GRCh38_alldifficultregions.bed.gz
giab_lowmap_segdup	$GIAB/Union/GRCh38_alllowmapandsegdupregions.bed.gz
giab_lowmappability	$GIAB/Mappability/GRCh38_lowmappabilityall.bed.gz
giab_segdups	$GIAB/SegmentalDuplications/GRCh38_segdups.bed.gz
giab_tandem_homopolymer	$GIAB/LowComplexity/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz
giab_satellites	$GIAB/LowComplexity/GRCh38_satellites_slop5.bed.gz
giab_gc_extreme	$GIAB/GCcontent/GRCh38_gclt25orgt65_slop50.bed.gz
giab_gaps	$GIAB/OtherDifficult/GRCh38_gaps_slop15kb.bed.gz
giab_other_difficult	$GIAB/OtherDifficult/GRCh38_allOtherDifficultregions.bed.gz
giab_false_dup_incorrect	$GIAB/OtherDifficult/GRCh38_false_duplications_incorrect_copy.bed.gz
giab_false_dup_correct	$GIAB/OtherDifficult/GRCh38_false_duplications_correct_copy.bed.gz
giab_collapsed_dup	$GIAB/OtherDifficult/GRCh38_collapsed_duplication_FP_regions.bed.gz
giab_population_cnv_fp	$GIAB/OtherDifficult/GRCh38_population_CNV_FP_regions.bed.gz
giab_l1h	$GIAB/OtherDifficult/GRCh38_L1H_gt500.bed.gz
giab_vdj	$GIAB/OtherDifficult/GRCh38_VDJ.bed.gz
giab_mhc	$GIAB/OtherDifficult/GRCh38_MHC.bed.gz
giab_kir	$GIAB/OtherDifficult/GRCh38_KIR.bed.gz
encode_blacklist	https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz
ucsc_rmsk	$UCSC/rmsk.txt.gz
ucsc_centromeres	$UCSC/centromeres.txt.gz
EOF
)

mkdir -p "$TRACK_DIR/giab_md5"
GIAB_MD5=$TRACK_DIR/giab_md5/GRCh38-genome-stratifications-md5s.txt
if [ ! -s "$GIAB_MD5" ]; then
    curl -sS -L --fail -o "$GIAB_MD5" "$GIAB/GRCh38-genome-stratifications-md5s.txt"
fi

TABLE=$TRACK_DIR/tracks.tsv
printf 'name\tfile\turl\tdownloaded\tmd5\tmd5_check\n' > "$TABLE"

while IFS=$'\t' read -r NAME URL; do
    FILE=$(basename "$URL")
    DIR=$TRACK_DIR/$NAME
    mkdir -p "$DIR"
    if [ ! -s "$DIR/$FILE" ]; then
        echo "fetch $NAME: $URL"
        curl -sS -L --fail -o "$DIR/$FILE.part" "$URL"
        mv "$DIR/$FILE.part" "$DIR/$FILE"
    else
        echo "keep  $NAME: $DIR/$FILE"
    fi
    MD5=$(md5sum "$DIR/$FILE" | cut -d ' ' -f 1)
    CHECK=none
    if [[ $NAME == giab_* ]]; then
        EXPECTED=$(grep -F "/$FILE" "$GIAB_MD5" | head -n 1 | cut -d ' ' -f 1 || true)
        if [ -z "$EXPECTED" ]; then
            CHECK=not_listed
        elif [ "$EXPECTED" = "$MD5" ]; then
            CHECK=ok
        else
            CHECK=mismatch
        fi
    fi
    DATE=$(date -r "$DIR/$FILE" +%Y-%m-%d)
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$NAME" "$NAME/$FILE" "$URL" "$DATE" "$MD5" "$CHECK" >> "$TABLE"
done <<< "$TRACKS"

echo "Tracks in $TRACK_DIR (list: $TABLE)"
if grep -q -P '\tmismatch$' "$TABLE"; then
    echo "md5 mismatch for some GIAB files; see $TABLE" >&2
    exit 1
fi
