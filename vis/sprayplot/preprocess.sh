#!/bin/bash

INFILE=output
OUTFILE=points

cat "$INFILE" | \
    sed -n '/^@@ __EVENT_SOURCE_TEST__$/ { s///; :a; n; p; ba; }' \
    > $OUTFILE
#   ^^ Returns all lines AFTER the matched pattern
#      s/// clears pattern space (apparently)
#      :a is a lable
#      n grabs next line
#      p prints it
#      ba says go back to :a
#

## Now output a copy that negates the longitude, effectively a
## reflection about the lon=0|180 plane, or an inversion of the east
## and west directions (since the Aki reference frame has +X pointing
## north).  This is used to plot beachballs which are kindof an
## "inside out" view.  Also, strip all but the P amplitudes.

cat "$OUTFILE" | \
    grep 'c$' | \
    awk '{ {$1*=-1} print $1, $2, $3, $4 }' > \
    "$OUTFILE"_refl

