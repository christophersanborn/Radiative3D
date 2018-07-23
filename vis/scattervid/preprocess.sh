#!/bin/bash
INFILE="reports.dat"
OUT_P="scatters_P"
OUT_S="scatters_S"
OUT_SRC="source_loc"

COL_TIME=5
COL_X=9
COL_Y=10
COL_Z=11

if [ ! -f $INFILE ]; then
  echo "No input files to process."
  exit
fi

# Get just the scatter events resulting in P waves:
cat $INFILE | \
    grep "^\(SCT\|REF\):[0-9 ]*P" | \
    # And retain just the relevant columns:
    awk '{ print $'$COL_TIME', $'$COL_X', $'$COL_Y', $'$COL_Z' }' > \
    $OUT_P

# Get just the scatter events resulting in S waves:
cat $INFILE | \
    grep "^\(SCT\|REF\):[0-9 ]*S" | \
    # And retain just the relevant columns:
    awk '{ print $'$COL_TIME', $'$COL_X', $'$COL_Y', $'$COL_Z' }' > \
    $OUT_S

# Get location of event source:
cat $INFILE | \
    grep -m 1 "^GEN:" | \
    # And retain just the relevant columns:
    awk '{ print $'$COL_X', $'$COL_Y', $'$COL_Z' }' > \
    $OUT_SRC

# Retrieve header text from INFILE:  (Every line up to the first GEN: event)
#awk '/^GEN:/ {exit} {print}' $INFILE > data_header

# Compress INFILE for archival:
gzip $INFILE

# Make framecache directory if it does not already exist
if [ ! -e "framecache" ]; then
    mkdir framecache
fi
