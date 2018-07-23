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

#SOURCETYP=$(EventParams "Selby" 0)     # Lookup source parameters
SOURCETYP=SDR,22.5,90,0                 # Or specify directly
SOURCELOC=0,0,-10                       # Source location
                              # Source param choices include 'expl',
                              # 'eq', 'Selby', etc.  Second paramter
                              # is usually isofrac (range [-1.0,1.0]
                              # or [-90,90])

FREQ=2.0                      # Phonon frequency to model
NUMPHONS=60K                  # Number of phonons to spray
RECTIME=350                   # Recording duration of seismometers.
GATHER=40.0                   # Terminal gather radius, in kilometers.

SCAT1=0.8,0.01,0.20,0.2,200   # Scat Args (nu,eps,a,kappa) and Q in sedi's
SCAT2=0.8,0.01,0.20,0.3,1500  # Scat Args in crust
SCAT3=0.8,0.01,0.20,0.3,1500  # Scat Args in crust, pinched region
SCAT4=0.8,0.01,0.20,0.4,1500  # Scat Args in moho
SCAT5=0.8,0.01,0.20,0.5,900   # Scat Args in mantle layer
                              # (Note: Q's specified are Q_s values.
                              # Q_p is computed from assumption that
                              # Q_kappa is infinite.)
LAYERS=2.0,30.0,5.0           # SediThick,CrustThick,MohoThick
PINCH=.3666667,.4736842,1,1   # PinchFrac,DepthFrac,SediFrac,MohoFrac

COMPARGS=$SCAT1,$SCAT2,$SCAT3,$SCAT4,$SCAT5,$LAYERS,$PINCH

MFPOVERRIDE=--overridemfp=25,50     # Uncomment to pin P,S Mean-Free-Paths.
NODEFLECT=--nodeflect               # Uncomment to make scattering
                                    # non-delfectionary.

#FLATTEN="--flatten"          # If set to "--flatten", apply Earth-flattening
FLATTEN=""                    # transformation to depths and velocities.
                              # (Set to "" to disable.)
ADDITIONAL=""                 # Additional params. (Such as --ocsraw,
                              # or --earthrad=xxxx) Radii: Earth: 6371 km,
                              # Mars: 3389, Mercury: 2440, Earth's
                              # moon: 1737

SEISORIG1=0,67.5,0            # Array locations: (Range, Azimuth, Z)
SEISORIG2=0,112.5,0           #
SEISORIG3=0,90,0              #
SEISDEST1=950,67.5,0          #
SEISDEST2=950,112.5,0         #
SEISDEST3=950,90,0            #
SEIS1=--seis-p2p=$SEISORIG1,$SEISDEST1,1.0,2.0,$GATHER,16
SEIS2=--seis-p2p=$SEISORIG2,$SEISDEST2,1.0,2.0,$GATHER,16
SEIS3=--seis-p2p=$SEISORIG3,$SEISDEST3,1.0,2.0,$GATHER,16
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
  STA="STA"
  STB="STB"
  STC="STC"
  MODELCODE="PINCH"
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
  scattervid_p2p(0,
            struct(
              "destination", [610,0,0],
              "seispat", "seisfiles/seis_%03d.octv",
              "seisrange", [0:32],
              "loctext", "Norwich, UK",
              "numframes", 300
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
              "seisrange", [0:48],
              "loctext", "Norwich, UK",
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
