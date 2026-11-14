#!/bin/bash
#
# Prism SETUP
#
# Run once after cloning:
#   ./Prism/00.SETUP.sh
#

set -eu

PDIR=$(dirname "$(readlink -f "$0")")

#
# Install nextflow in Prism/bin
#
# conf/ENV.sh sets NXF_VER (the pinned version) and NXF_HOME, so the
# install does not write to the home folder.
#

source "$PDIR/conf/ENV.sh"

# Never write to the node's /tmp
mkdir -p /localscratch/soccin/oncoanalyser
SETUP_TMP=$(mktemp -d -p /localscratch/soccin/oncoanalyser setup_XXXXXX)
export TMPDIR=$SETUP_TMP
export NXF_OPTS="$NXF_OPTS -Djava.io.tmpdir=$SETUP_TMP"

mkdir -p "$PDIR/bin"
cd "$PDIR/bin"
curl -s https://get.nextflow.io | bash

echo
echo "Installed nextflow $NXF_VER in $PDIR/bin"
