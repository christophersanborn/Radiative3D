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
##  This DO-SCRIPT sets up a run based on the Crust Pinch WCG model
##

## One-liner description: (Keep this BRIEF.)
##
INTENT="Example crust pinch simulation."
CAMPAIGN="" ## (e.g. "For AGU Poster Dec 2016" or some such.)

SIMTARGET="waveform"          # Choice: 'waveform' or 'video'. Affects
                              # defaults not otherwise specified.
MODIDX=5  # NSCP              # Model Index: Selects from custom coded models.
                              # 1: Lop Nor (base or moho depends on COMPARGS),
                              # 5: North Sea Crust Pinch model,
                              # 8: Crust Upthrust model

#SOURCETYP=$(EventParams "Selby" 0)     # Lookup source parameters
SOURCETYP=SDR,22.5,90,0                 # Or specify directly
SOURCELOC=0,0,-10                       # Source location
                              # Source param choices include 'expl', 'eq',
                              # 'Selby', etc.  Second paramter is usually
                              # isofrac (if range [-1.0,1.0]) or isoangle
                              # (if range [-90,90]).

FREQ=2.0                      # Phonon frequency to model
NUMPHONS=1M                   # Number of phonons to spray (Recommend: 50M)
RECTIME=600                   # Recording duration of seismometers (seconds).
BINSIZE=2.00                  # Seismometer time-bin size in seconds
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
SEIS1=--seis-p2p=$SEISORIG1,$SEISDEST1,1.0,2.0,$GATHER,160
SEIS2=--seis-p2p=$SEISORIG2,$SEISDEST2,1.0,2.0,$GATHER,160
SEIS3=--seis-p2p=$SEISORIG3,$SEISDEST3,1.0,2.0,$GATHER,160


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
## GENERATE FIGURES AND GRAPHICS:       #   be carved off into 'do-figsonly.sh'
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
if [ -x "`which pdfcrop 2> /dev/null`" ]; then  # Check whether pdfcrop is
pdfcroptobb() {  # $1=tmp-file $2=ouputfile     # installed. If yes, then 
  pdfcrop --margins 0 "$1" "$2" >/dev/null  # define pdfcroptobb (PDF-Crop-  
  rm "$1"                                   # to-Bounding-Box), otherwise
} else                                      # define a dummy function.
pdfcroptobb() {                             #
  echo "pdfcrop not installed."; 
  mv "$1" "$2"; }
fi


# Make seismometer and model maps:
echo "Plotting model maps and cross-sections."
cat stdout.txt | \
    grep -A 10000 "#  R3D_GRID:" | \
    grep -B 10000 "#  END R3D_GRID" | \
    revtail -n +2 | \
    tail -n +9 | \
    sed 's/\*\*\*/000/g' > griddump.txt  
octave -qf <<"EOF"
  modelplot_WCGrangeelev(
      "griddump.txt", 0,
      "Crust Profile with Pinch Region",
      "out_mparams.octv",
      "seisfiles/seis_%03d.octv",
      [0,9:10:159]);
  axis([0 950]);
  print("profile-model.png","-r200");   # PNG version for wiki
  print("tmp-profile-model.pdf");       # PDF version for publication
EOF
pdfcroptobb "tmp-profile-model.pdf" "profile-model.pdf"

# Seismograms:
if [ "$1" != "noseis" ]; then
echo "Plotting envelopes..."
for file in seisfiles/{seis_000.octv,seis_?[13579][9].octv} ; do
    ofile=`basename "${file%.octv}.png"`  # PNG for wiki
    echo $ofile
    octave -qf << EOF
      seisplot("$file"); 
      papertext(0.01, 0.00, "$RUNID",   ## Plot RunID in lowerleft corner
                "verticalalignment", "bottom", "horizontalalignment", "left",
                "fontsize", 5, "fontweight", "normal", "color", [0 0 0]);
      print("$ofile","-r200");
EOF
done
fi
for file in seisfiles/{seis_159.octv,seis_319.octv,seis_479.octv} ; do
    ofile=`basename "${file%.octv}.pdf"`  # PDF for publication
    echo $ofile
    octave -qf << EOF
      seisplot("$file"); 
      hRID = papertext(0.01, 0.00, "$RUNID",   ## Plot RunID in lowerleft corner
                  "verticalalignment", "bottom", "horizontalalignment", "left",
                  "fontsize", 5, "fontweight", "normal", "color", [1 1 1]);
      #set(hRID, "color", [1 1 1]);  # Invisible but selectable in pdf version
      #title("");  # Kill title, and Kill 1/4-inch extra margin it claims:
      #set(gcf(),"paperposition", get(gcf(),"paperposition")-[0 0 0 0.25])
      print("tmp-$ofile");
EOF
    pdfcroptobb "tmp-$ofile" "$ofile"  # Crop to bb
done

# Traveltime Curves (Template Function):
produce_ttcurves() {  # $1:     station code
                      # $2, $3: seis indices
                      # $4:     caxis limit
  # Create links to reference curves. User must manually define dir-
  # ectory links "other..." to activate.  Facilitates energy compar-
  # isons to a baseline run. (Link "other" to the baseline run,
  # "other2" and "other3" to other test-case runs.)
  ncbase="traveltime-$1-xyz-normcurve"
  [ -d other ] && ln -sf "other/${ncbase}_smooth.octv" "${ncbase}_common.octv"
  [ -d other ] && ln -sf "other/${ncbase}.octv" "${ncbase}_other.octv"
  [ -d other2 ] && ln -sf "other2/${ncbase}.octv" "${ncbase}_other2.octv"
  [ -d other3 ] && ln -sf "other3/${ncbase}.octv" "${ncbase}_other3.octv"
  octave -qf <<EOF
    AR = array("seisfiles/seis_%03d.octv",$2,$3);
    NORMCURVE = arrayimage(AR, [1 1 1], 2.0, [0 400]);
    save("traveltime-$1-xyz-normcurve.octv", "NORMCURVE");
    NC2 = NORMCURVE;
    NORMCURVE = normcurve_fitpowerlaw(NORMCURVE,10,160);
    save("traveltime-$1-xyz-normcurve_smooth.octv", "NORMCURVE");
    if (exist("traveltime-$1-xyz-normcurve_common.octv","file")==2)
      NORMCURVE = load("traveltime-$1-xyz-normcurve_common.octv").NORMCURVE;
      if (exist("traveltime-$1-xyz-normcurve_other.octv","file")==2)
        NC1 = load("traveltime-$1-xyz-normcurve_other.octv").NORMCURVE;
        normcurve_compare(NC1, NC2, NORMCURVE);
        print("traveltime-$1-xyz-etracecompare.png", "-r160");
      end 
    end   # Loads a common curve if one has been created (requires
          # user intervention).  If _othercurve also exists, 
          # then plot a comparisson.
    ROI{1}.RegionWindow = [310 530];    # Region of Interest definitions
    ROI{1}.RegionVSpan = [0.04 0.51];   #
    ROI{1}.RegionLineWidth = 1.0;
    ROI{1}.Label = "";
    ROI{1}.LabelFontSize = 8;
    ROI{1}.LabelVPos = 0.10;
    ROI{2}.RegionWindow = [370 470];
    ROI{2}.RegionVSpan = [0.09 0.46];
    ROI{2}.RegionLineWidth = 1.0;
    ROI{2}.Label = "";
    ROI{2}.LabelFontSize = 8;
    ROI{2}.LabelVPos = 0.35;
    arrayimage(AR, [1 1 1], 2.0, [0 400], NORMCURVE, $4);
    annotate_regions(ROI);
    annotate_array(); # (mark phases with velocity lines)
    hRID = papertext(0.03, 0.00, "$RUNID",   ## Plot RunID in lowerleft corner
              "verticalalignment", "bottom", "horizontalalignment", "left",
              "fontsize", 5, "fontweight", "normal");
    #set(findobj("tag", "colorbar"), "yticklabel", []); # kill colorbar labels
    #set(findobj("tag", "colorbar"), "ytick", []);      # kill clrbar ticks too
    print("traveltime-$1-xyz-g2.0.png", "-r200");       ## PNG version **
    #set(hRID, "color", [1 1 1]);  # Invisible but selectable in pdf version
    #title("");  # Kill title
    print("tmp-traveltime-$1-xyz-g2.0.pdf");            ## PDF version **
EOF
    pdfcroptobb "tmp-traveltime-$1-xyz-g2.0.pdf" \
                "traveltime-$1-xyz-g2.0.pdf"  # Crop to bb
}
# TT Curves:
produce_ttcurves STA   0 159 8    # (Out the S-wave principle axis)
produce_ttcurves STB 160 319 8    # (Out the P-wave principle axis)
produce_ttcurves STC 320 479 8    # (Out the equivalency axis)

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
