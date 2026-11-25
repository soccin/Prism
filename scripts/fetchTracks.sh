#!/bin/bash
# Prism: download the GRCh38 region tracks that postOA.sh uses to
# annotate SV breakends with technical context (mappability, repeats,
# segmental duplications, blacklists, gaps). Run once per site, on the
# login node (downloads only):
#   ./Prism/scripts/fetchTracks.sh [TRACK_DIR]
# TRACK_DIR defaults to /scratch/core001/bic/socci/reference/prism_tracks/GRCh38.
# Each track goes into its own folder TRACK_DIR/<name>/. TRACK_DIR/tracks.tsv
# lists every track: name, file, source URL, download date, md5, and the
# result of the md5 check against the md5 the source publishes (GIAB
# list, gnomAD bucket, 1000 Genomes manifest); none where the source
# publishes none. Files already present are kept and not downloaded
# again; delete a track's folder to fetch it again.
# The population SV VCFs (gnomad_sv, kgp_sv, dbvar_*) are then reduced to
# the fields svCohort.R uses by scripts/popSvSites.sh, a job (see that
# script).

set -eu
set -o pipefail

TRACK_DIR=${1:-/scratch/core001/bic/socci/reference/prism_tracks/GRCh38}
TRACK_DIR=${TRACK_DIR%/}

GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all
UCSC=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database
# gnomAD-SV v4.1 sites (Google Cloud public bucket). The md5 is checked
# against the one the bucket lists for the object.
GNOMAD_BUCKET=gcp-public-data--gnomad
GNOMAD_SV=release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz
GNOMAD=https://storage.googleapis.com/$GNOMAD_BUCKET
# Population SV call sets from short-read WGS only (no long-read sets):
#   kgp_sv         1000 Genomes high coverage, 3,202 samples (Byrska-Bishop
#                  2022); md5 from the release manifest
#   dbvar_nstd166  gnomAD-SV v2.1, 14,891 genomes (Collins 2020), as dbVar
#                  study nstd166; dbVar lists no md5
#   dbvar_nstd200  CCDG, 17,795 genomes (Abel 2020), as dbVar study
#                  nstd200; dbVar lists no md5
KGP=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20210124.SV_Illumina_Integration
KGP_MD5=1KGP_3202.Illumina_ensemble_callset.freeze_V3.manifest.md5
DBVAR=ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf

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
gnomad_sv	$GNOMAD/$GNOMAD_SV
kgp_sv	$KGP/1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz
dbvar_nstd166	$DBVAR/nstd166.GRCh38.variant_call.vcf.gz
dbvar_nstd200	$DBVAR/nstd200.GRCh38.variant_call.vcf.gz
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
    # VCF tracks: also fetch the tabix index (popSvSites.sh checks the
    # record count against it).
    if [[ $FILE == *.vcf.gz ]] && [ ! -s "$DIR/$FILE.tbi" ]; then
        echo "fetch $NAME: $URL.tbi"
        curl -sS -L --fail -o "$DIR/$FILE.tbi.part" "$URL.tbi"
        mv "$DIR/$FILE.tbi.part" "$DIR/$FILE.tbi"
    fi
    MD5=$(md5sum "$DIR/$FILE" | cut -d ' ' -f 1)
    # The md5 the source publishes; empty: the source lists none.
    HAS_LIST=yes
    case $NAME in
        gnomad_sv)
            OBJECT=${GNOMAD_SV//\//%2F}
            EXPECTED=$(curl -sS -L --fail "https://storage.googleapis.com/storage/v1/b/$GNOMAD_BUCKET/o/$OBJECT?fields=md5Hash" \
                | sed -n 's/.*"md5Hash": *"\([^"]*\)".*/\1/p' | base64 -d | od -An -tx1 | tr -d ' \n')
            ;;
        kgp_sv)
            # The manifest has CRLF line ends.
            EXPECTED=$(curl -sS -L --fail "$KGP/$KGP_MD5" | tr -d '\r' | awk -F'\t' -v f="$FILE" '$1 == f {print $3}')
            ;;
        giab_*)
            EXPECTED=$(grep -F "/$FILE" "$GIAB_MD5" | head -n 1 | cut -d ' ' -f 1 || true)
            ;;
        *)
            HAS_LIST=no
            EXPECTED=""
            ;;
    esac
    if [ "$HAS_LIST" = no ]; then
        CHECK=none
    elif [ -z "$EXPECTED" ]; then
        CHECK=not_listed
    elif [ "$EXPECTED" = "$MD5" ]; then
        CHECK=ok
    else
        CHECK=mismatch
    fi
    DATE=$(date -r "$DIR/$FILE" +%Y-%m-%d)
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$NAME" "$NAME/$FILE" "$URL" "$DATE" "$MD5" "$CHECK" >> "$TABLE"
done <<< "$TRACKS"

echo "Tracks in $TRACK_DIR (list: $TABLE)"
if grep -q -P '\tmismatch$' "$TABLE"; then
    echo "md5 mismatch for some files; see $TABLE" >&2
    exit 1
fi
