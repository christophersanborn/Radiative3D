#!/bin/bash
##
##  "DO"-script for a Radiative3D run:
##
##  Copy and edit this file to easily specify (and remember)
##  parameters used for a particular run.
##
##  If $1 is set, it is used as an extra identifier token in the
##  output directory name.
##

##
##  This DO-SCRIPT sets up a run based on the LOP NOR cylinder model
##

R3D_EXE=../../main  # Executable to use
OUT_BASE=data       # Parent directory in which to put output
                    # subdirectories

event=eq        # See case statement below for choices


##
##  Lop Nor geography and event source parameters:
##

LOPXYZ=492.31,-263.65,1.05    # Lop Nor, Surface, site of Chinese Test 596
XINXYZ=425.54,-169.53,0.98    # Lop Nor, Surface, site of Xinjiang Quake 030313
MAKXYZ=-102.27,430.84,0.60    # Station MAK
WUSXYZ=-390.04,-167.18,1.457  # Station WUS

FREQ=2.0                      # Phonon frequency to model
NUMPHONS=200K                 # Number of phonons to spray
RECTIME=600                   # Recording duration of seismometers.
CYLRAD=1200                   # Total (cylindrical) radius of model.
                              # (Phonons exceeding this radius from XY
                              # = (0,0), or this travel time from t=0,
                              # are abandoned.)

SCAT1=0.5,0.01,0.5,0.2        # Scattering Args (nu,eps,a,kappa) in sediments
SCAT2=0.5,0.01,0.5,0.3        # Scattering Args (nu,eps,a,kappa) in crust
SCAT3=0.5,0.01,0.7,0.5        # Scattering Args (nu,eps,a,kappa) in mantle


#
#   Event Parameters:
#

case "$event" in
    expl)   # Generic explosion
        SOURCELOC=425.54,-169.53,-1.02    # Lop Nor, 2km below surface
        SOURCETYP=EXPL                    #
        ;;
    eq)     # Xinjiang earthquake
        SOURCELOC=425.54,-169.53,-31.02   # Xinjiang, 32km below surface
        SOURCETYP=SDR,125,40,90           #
        ;;
    *)
        echo Event code not recognized.
        exit
        ;;
esac


## PREPROCESS INPUT FILES:
##
## (If any pre-processing is needed, code it here.)
##

#if ! make -q
#then
#    echo Make indicates $R3D_EXE is not up to date.
#    echo Please run make and re-run this script.
#    exit
#fi

#
#  Set up directory in which to put output and input records:
#
#           Directory will be named with a timestamp and an
#           extra token if $1 is set.
#

tstamp=`date +"%Y%m%d-%H%M%S"`
outdirname="$tstamp"-"$1"-R3D
outdir="$OUT_BASE/$outdirname"
mkdir -p "$outdir"

if [ -L "$OUT_BASE/previous" ]; then    # Generate "latest" and "previous"
    rm "$OUT_BASE/previous"             # links to make it easier to cd
fi                                      # into new output
if [ -L "$OUT_BASE/latest" ]; then
    mv "$OUT_BASE/latest" "$OUT_BASE/previous"
fi
ln -sf "$outdirname" "$OUT_BASE/latest"
#

cp "$0" "$outdir"/                # Copy this file
cp $R3D_EXE "$outdir"/            # Keep the exe for replicability
R3D_EXE="$outdir"/$(basename $R3D_EXE)  # Run the copy instead of the orig
#svn info > "$outdir"/svn-info.txt # Keep revision info
#svn st >> "$outdir"/svn-info.txt  #
#svn diff >> "$outdir"/svn-info.txt 
cp *.m preprocess.sh "$outdir"/    # We will use these in post-process


LOGFILE="$outdir"/logfile
> $LOGFILE


## RUN THE SIMULATION:
##
## (Edit command-line args here.)
##

echo Begin Radiative3D Run: `date +"%Y.%m.%d-%H:%M:%S"` >> "$LOGFILE"

$R3D_EXE --reports=ALL_ON \
        --output-dir="$outdir" \
        --report-file=reports.dat \
        --mparams-outfile=out_mparams.octv \
        --dump-grid \
        --frequency=$FREQ \
        --num-phonons=$NUMPHONS \
        --toa-degree=8 \
        --source=$SOURCETYP \
        --source-loc=$SOURCELOC \
        --timetolive=$RECTIME \
    | tee "$outdir"/stdout.txt

echo End__ Radiative3D Run: `date +"%Y.%m.%d-%H:%M:%S"` >> "$LOGFILE"

## POST-PROCESS OUTPUT FILES:
##
## (If any post-processing is needed, code it here.)
##

echo Begin Figure Generation: >> "$LOGFILE"

#
# Scatter Videos:
#

rwd=`pwd` 
cd "$outdir"
echo "Processing output for visualization:"
./preprocess.sh

# Vid from above
echo "Making Scatter-Vid Above View:"
octave -qf <<EOF
  scattervid_above();
EOF
mv framecache/movie.mp4 scatvid-above.mp4
mv framecache framecache-above
mkdir framecache

# Vid from elevation viewpoint
echo "Making Scatter-Vid Elevation View:"
octave -qf <<EOF
  #scattervid_p2p();
  scattervid_p2p([1000,1000,0]);  # Specify dest for auto-range
EOF
mv framecache/movie.mp4 scatvid-elev.mp4
mv framecache framecache-elev
mkdir framecache

cd "$rwd"

echo End Figure Generation: `date +"%Y.%m.%d-%H:%M:%S"` >> "$LOGFILE"

## END
##
