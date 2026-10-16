#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
threads=$1
outdir=$2
phonons=${3:-24}
mkdir -p "$outdir"

# This is the simulation-only portion of do-lopnor.sh's Lop Nor waveform case.
# The production do-script continues into Octave figure generation after this
# command; this test intentionally stops after the executable returns.
LOPXYZ=492.31,-263.65,1.05
XINXYZ=425.54,-169.53,0.98
MAKXYZ=-102.27,430.84,0.60
WUSXYZ=-390.04,-167.18,1.457
SOURCELOC=425.54,-169.53,-31.02
SOURCETYP=SDR,125,40,90,0.0
SCAT1=0.8,0.01,0.5,0.2,50
SCAT2=0.8,0.01,0.5,0.3,1000
SCAT3=0.8,0.01,0.7,0.5,300
COMPARGS=$SCAT1,$SCAT2,$SCAT3

set -- "$ROOT/main" \
  --reports=ALL_OFF \
  --output-dir="$outdir" \
  --report-file=reports.dat \
  --mparams-outfile=out_mparams.octv \
  --num-phonons="$phonons" \
  --toa-degree=5 \
  --source="$SOURCETYP" \
  --source-loc="$SOURCELOC" \
  --frequency=2.0 \
  --timetolive=0.02 \
  --binsize=0.001 \
  --grid-compiled=1 \
  --range=1200 \
  --flatten \
  --model-args="$COMPARGS" \
  --seed=596 \
  --seis-p2p="$XINXYZ,$WUSXYZ,1.0,2.0,40.0,8" \
  --seis-p2p="$XINXYZ,$MAKXYZ,1.0,2.0,40.0,8"

if [ "$threads" != auto ]; then
  set -- "$@" --threads="$threads"
fi

(
  cd "$outdir"
  "$@" >stdout.txt 2>stderr.txt
)

# Guard the test contract: the comparison runner must never invoke Octave.
if grep -i 'octave' "$outdir/stdout.txt" "$outdir/stderr.txt" >/dev/null 2>&1; then
  echo "unexpected Octave invocation/output in simulation-only Lop Nor case" >&2
  exit 1
fi
