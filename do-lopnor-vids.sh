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
##  (Note: Do not edit line 3.  It is used by subsequent scripts to
##  identify this file as a "do script". (See wikify.sh.))
##
source scripts/do-fundamentals.sh
##
##  This DO-SCRIPT sets up a VIDEO run based on the Crust Pinch WCG model
##

SIMTARGET="video"             # Choice: 'waveform' or 'video'. Affects
                              # defaults not otherwise specified.
MODIDX=21  # LOP NOR           # Model Index: Selects from custom coded models.
                              # 1: Lop Nor (base or moho depends on COMPARGS),
                              # 5: North Sea Crust Pinch model,
                              # 8: Crust Upthrust model

event=eq        # 'eq' or 'expl' - See case statement below for choices
EQISOFRAC=0.0   # Iso fraction for EQ event (choose from range [-1.0, 1.0])

##
##  Lop Nor geography and event source parameters:
##

LOPXYZ=492.31,-263.65,1.05    # Lop Nor, Surface, site of Chinese Test 596
XINXYZ=425.54,-169.53,0.98    # Lop Nor, Surface, site of Xinjiang Quake 030313
MAKXYZ=-102.27,430.84,0.60    # Station MAK
WUSXYZ=-390.04,-167.18,1.457  # Station WUS

FREQ=2.0                      # Phonon frequency to model
NUMPHONS=3K                  # Number of phonons to spray
RECTIME=350                   # Recording duration of seismometers.
GATHER=40.0                   # Terminal gather radius, in kilometers.
CYLRAD=1200                   # Total (cylindrical) radius of model.
                              # (Phonons exceeding this radius from XY
                              # = (0,0), or this travel time from t=0,
                              # are abandoned.)

SCAT1=0.8,0.01,0.5,0.2,50     # Scat Args (nu,eps,a,kappa) and Q in sedi's
SCAT2=0.8,0.01,0.5,0.3,1000   # Scat Args (nu,eps,a,kappa) and Q in crust
SCAT3=0.8,0.01,0.7,0.5,300    # Scat Args (nu,eps,a,kappa) and Q in mantle
                              # (Note: Q's specified are Q_s values.
                              # Q_p is computed from assumption that
                              # Q_kappa is infinite.)
# To enable Moho structure, uncomment the following:
SCAT3=$SCAT3,0.8,0.02,0.5,0.5,2000  # Scat Args in transition layers

COMPARGS=$SCAT1,$SCAT2,$SCAT3

MFPOVERRIDE=--overridemfp=1,1     # Uncomment to pin P,S Mean-Free-Paths.
NODEFLECT=--nodeflect               # Uncomment to make scattering
                                    # non-delfectionary.

FLATTEN="--flatten"          # If set to "--flatten", apply Earth-flattening
#FLATTEN=""                    # transformation to depths and velocities.
                              # (Set to "" to disable.)
ADDITIONAL=""                 # Additional params. (Such as --ocsraw,
                              # or --earthrad=xxxx) Radii: Earth: 6371 km,
                              # Mars: 3389, Mercury: 2440, Earth's
                              # moon: 1737

#
#   Event Parameters:
#

case "$event" in
    expl)   # Generic explosion
        SOURCELOC=425.54,-169.53,-1.02    # Lop Nor, 2km below surface
        #SOURCELOC=425.54,-169.53,-5.02    # Lop Nor, 6km below surface
        #SOURCELOC=425.54,-169.53,-31.02   # Xinjiang, 32km below surface
        SOURCETYP=EXPL                    #
        ;;
    eq)     # Xinjiang earthquake
        #SOURCELOC=425.54,-169.53,-1.02    # Lop Nor, 2km below surface
        #SOURCELOC=425.54,-169.53,-5.02    # Xinjiang, 6km below surface
        SOURCELOC=425.54,-169.53,-31.02   # Xinjiang, 32km below surface
        SOURCETYP=SDR,125,40,90,$EQISOFRAC  #
        ;;
    *)
        echo Event code not recognized.
        exit
        ;;
esac

SEISORIG1=$XINXYZ             # Array locations:
SEISORIG2=$XINXYZ             #
#SEISORIG3=$XINXYZ            #
SEISDEST1=$WUSXYZ             #
SEISDEST2=$MAKXYZ             #
#SEISDEST3=""                 #
SEIS1=--seis-p2p=$SEISORIG1,$SEISDEST1,1.0,2.0,$GATHER,16
SEIS2=--seis-p2p=$SEISORIG2,$SEISDEST2,1.0,2.0,$GATHER,16
#SEIS3=--seis-p2p=$SEISORIG3,$SEISDEST3,1.0,2.0,$GATHER,160
BINSIZE=10.0                  # Large time-bins, not useful for trace analy-
                              # sis, but makes seis files smaller. (We only
                              # want them to plot a map as the video
                              # backdrop.)

## RUN THE SIMULATION:
##
##  Set up output directory and run the sim:
##

PopDefaults $SIMTARGET        ## Defined in do-fundamentals.sh
CheckBuildStatus              ##  ''
CreateOutputDirectory $@      ##  ''
RunSimulation                 ##  ''

echo Begin Figure Generation: >> "$LOGFILE"
rwd=`pwd`       # Switch to output directory - the rest of our work
cd "$outdir"    # will be done there.
## (Do not edit/remove this comment block.)
#!/bin/bash
## ___FIG_GEN_START___
##

##                                      # Note: this part of the do-script will
## GENERATE VIDEOS:                     #   be carved off into 'do-figsonly.sh'
##                                      #   to facilitate re-running and/or
                                        #   customizing figure output after the
                                        #   initial run.

# Organize seisfiles
mkdir -p seisfiles  # Stash seis files in a subdirectory (aesthetics)
mv seis_???.octv seisfiles/ 2> /dev/null || : # (suppress error msg and code)

# Text values for Wiki script:
cat > wikify.incl <<EOF
  STA="WUS"
  STB="MAK"
  STC="xxx"
  MODELCODE="LOP"
EOF

# Determine RunID from name of do-script:   (for labeling figs)
SHFILE=RUNID
for file in *.sh; do
    if grep -qx '##  "DO"-script for a Radiative3D run:' "$file"
    then                # Catch ONLY file with R3D do marker
        SHFILE="$file"  # Also catches no more than one do- file
    fi                  #
done
RUNID=${SHFILE%.sh}
RUNID=${RUNID#do-}
echo "Found RunID: $RUNID"

# Some octave/gnuplot set-up and platform-independence stuff:
alias octave="octave --no-window-system"    # Prevents trouble on some systems.
export GNUTERM="dumb"                   # Ascii, but we don't disply so is OK.
shopt -s expand_aliases                 # Aliases in scripts... needs enabled.
reverser="tac"                                  # Linux has tac (backwords cat)
[ "`uname`" == "Darwin" ] && reverser="tail -r" # But OS X has 'tail -r',...
revtail() { $reverser | tail $@ | $reverser; }  # 'head' with neg line count
                                                # does not work on OS X, so
                                                # use revtail in its place.

##
## Pre-process:
##
echo "Processing output for visualization:"
./preprocess.sh
cat stdout.txt | \
    grep -A 10000 "#  R3D_GRID:" | \
    grep -B 10000 "#  END R3D_GRID" | \
    revtail -n +2 | \
    tail -n +9 | \
    sed 's/\*\*\*/000/g' > griddump.txt

##
## Elevation (Side-View) Video:
##
echo "Making Scatter-Vid Elevation View:"
mkdir -p framecache
octave -qf <<EOF
  scattervid_model(0,
            struct(
              "destination","seisfiles/seis_031.octv", # 
              "gridfile",  "griddump.txt",             # A 3x1xN Cylinder grid
              "gridpair",  [1 2],                      # (1=LOP, 2=MAK, 3=WUS)
              "mparams",   "out_mparams.octv",         # Metadata
              "window",    [-0.05 1.4 10 -130],        # [LRTB], L,R fractional
              "seispat",   "seisfiles/seis_%03d.octv", # Seismometers files
              "seisrange", [],                         # ...which ones to plot
              "txtkey",    {{"LOP","MAK","WUS"}},      # Location names
              "loctext",   "Lop Nor",                  # Source verbose name
              "base",      [492.310, -263.650, 1.050], # Location of LOP
              "duration",  600,                        # Time, (if<PhononTTL)
              "numframes", 300                         # Total frames in video
              )
            );
EOF
mv framecache/movie.mp4 scatvid-elev.mp4
if [ -d framecache-elev ]; then
    rm -rf framecache-elev-old
    mv framecache-elev framecache-elev-old
fi
mv framecache framecache-elev

##
## Plan-View (from above) Video:
##
echo "Making Scatter-Vid Above View:"
mkdir -p framecache
octave -qf <<EOF
  scattervid_above(0,
            struct(
              "seispat",   "seisfiles/seis_%03d.octv",
              "seisrange", [0:32],
              "loctext",   "Lop Nor",
              "numframes", 300
            )
          );
EOF
mv framecache/movie.mp4 scatvid-above.mp4
if [ -d framecache-above ]; then
    rm -rf framecache-above-old
    mv framecache-above framecache-above-old
fi
mv framecache framecache-above


## (Do not edit/remove this comment block.)
## ___FIG_GEN_END___
##
cd "$rwd"
echo End Figure Generation: `date +"%Y.%m.%d-%H:%M:%S"` >> "$LOGFILE"
echo '///'
echo "Figure generation complete. To re-run figures, run 'do-figsonly.sh' in"
echo "the output directory, which may be modified to customize figure output."
echo "Output has been placed in $outdir."
echo "Logfile contents:"
cat "$outdir"/logfile
## END
##
