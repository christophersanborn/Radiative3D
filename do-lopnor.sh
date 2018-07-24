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
##  This DO-SCRIPT sets up a run based on the LOP NOR cylinder model
##

## One-liner description: (Keep this BRIEF.)
##
INTENT="Example Lop Nor simulation."
CAMPAIGN="" ## (e.g. "For AGU Poster Dec 2016" or some such.)

SIMTARGET="waveform"          # Choice: 'waveform' or 'video'. Affects
                              # defaults not otherwise specified.
MODIDX=1  # LOP NOR           # Model Index: Selects from custom coded models.
                              # 1: Lop Nor (base or moho depends on COMPARGS),
                              # 5: North Sea Crust Pinch model,
                              # 8: Crust Upthrust model
                              # 40: Halfspace model

event=eq        # 'eq' or 'expl' - See case statement below for choices
EQISOFRAC=0.0   # Iso fraction for EQ event (choose from range [-1.0, 1.0])
                # (Or choose as Iso Angle in range [-90, -1), (1, 90].)

##
##  Lop Nor geography and event source parameters:
##

LOPXYZ=492.31,-263.65,1.05    # Lop Nor, Surface, site of Chinese Test 596
XINXYZ=425.54,-169.53,0.98    # Lop Nor, Surface, site of Xinjiang Quake 030313
MAKXYZ=-102.27,430.84,0.60    # Station MAK
WUSXYZ=-390.04,-167.18,1.457  # Station WUS

FREQ=2.0                      # Phonon frequency to model
NUMPHONS=100K                 # Number of phonons to spray (Recommend: 140M)
RECTIME=600                   # Recording duration of seismometers (seconds).
BINSIZE=2.00                  # Seismometer time-bin size in seconds
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
#SCAT3=$SCAT3,0.8,0.02,0.5,0.5,2000  # Scat Args in transition layers

COMPARGS=$SCAT1,$SCAT2,$SCAT3

FLATTEN="--flatten"           # If set to "--flatten", apply Earth-flattening
#FLATTEN=""                   # transformation to depths and velocities.
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
SEIS1=--seis-p2p=$SEISORIG1,$SEISDEST1,1.0,2.0,$GATHER,160
SEIS2=--seis-p2p=$SEISORIG2,$SEISDEST2,1.0,2.0,$GATHER,160
#SEIS3=--seis-p2p=$SEISORIG3,$SEISDEST3,1.0,2.0,$GATHER,160


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
  seismap("seisfiles/seis_??9.octv");
  legend("hide"); # (We're not using separate P and S gather anymore)
  print("seismap.png","-r200");   # PNG version for wiki
  print("tmp-seismap.pdf");       # PDF version for publication
  ## LOP to MAK:
  modelplot("griddump.txt", [1 2], {"LOP","MAK"},
            "out_mparams.octv",
            "seisfiles/seis_%03d.octv",
            [160,169:10:319]);
  print("profile-LOP-MAK.png","-r200"); # PNG version for wiki
  #delete(findall("tag","sourcemark")); # Remove source marker for publication
  #delete(findall("tag","seismometer"));# Remove seismometer marks if desired
  print("tmp-profile-LOP-MAK.pdf");     # PDF version for publication
  ## LOP to WUS:
  modelplot("griddump.txt", [1 3], {"LOP","WUS"},
            "out_mparams.octv",
            "seisfiles/seis_%03d.octv",
            [0,9:10:159]);
  print("profile-LOP-WUS.png","-r200");
  ## WUS to MAK:
  modelplot("griddump.txt", [3 2], {"WUS","MAK"},
            "out_mparams.octv",
            "",
            []);
  print("profile-WUS-MAK.png","-r200");
EOF
pdfcroptobb "tmp-seismap.pdf" "seismap.pdf"  # Crop to bb
pdfcroptobb "tmp-profile-LOP-MAK.pdf" "profile-LOP-MAK.pdf"

# Seismograms:
if [ "$1" != "noseis" ]; then
echo "Plotting envelopes..."
for file in seisfiles/{seis_000.octv,seis_?[13579][9].octv} ; do
    ofile=`basename "${file%.octv}.png"`  # PNG's for wiki
    echo $ofile
    octave -qf << EOF
      seisplot("$file"); 
      papertext(0.01, 0.00, "$RUNID",   ## Plot RunID in lowerleft corner
                "verticalalignment", "bottom", "horizontalalignment", "left",
                "fontsize", 5, "fontweight", "demi", "color", [0 0 0]);
      print("$ofile","-r200");
EOF
done
fi
for file in seisfiles/{seis_319.octv,seis_159.octv} ; do
    ofile=`basename "${file%.octv}.pdf"`  # PDF for publication
    echo $ofile
    octave -qf << EOF
      seisplot("$file"); 
      hRID = papertext(0.01, 0.00, "$RUNID",   ## Plot RunID in lowerleft corner
                  "verticalalignment", "bottom", "horizontalalignment", "left",
                  "fontsize", 5, "fontweight", "demi", "color", [0 0 0]);
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
                      # $5:     caxis limit for 2nd plot (optional)
  # Create links to reference curves. User manually must define directory links
  # "other..." to activate. Enables energy comparisons to a baseline run. (Link
  # "other" to the baseline, "other2" and "other3" to other test-case runs.)
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
    ## VARIANT 1: Make one with norm relative to best-fit NORMCURVE: 
    arrayimage(AR, [1 1 1], 2.0, [0 400], NORMCURVE, $4);
    annotate_array();
    hRID = papertext(0.03, 0.00, "$RUNID",   ## Plot RunID in lowerleft corner
              "verticalalignment", "bottom", "horizontalalignment", "left",
              "fontsize", 5, "fontweight", "demi");
    print("traveltime-$1-xyz-g2.0.png", "-r200");
    ## VARIANT 2: And make one where norm relative to local peak and average:
    arrayimage(AR, [1 1 1], 2.0, [0 400], 0.5, ${5:-0});
    annotate_array(); # (mark phases with velocity lines)
    hRID = papertext(0.03, 0.00, "$RUNID",   ## Plot RunID in lowerleft corner
              "verticalalignment", "bottom", "horizontalalignment", "left",
              "fontsize", 5, "fontweight", "demi");
    set(findobj("tag", "colorbar"), "yticklabel", []); # kill colorbar labels
    set(findobj("tag", "colorbar"), "ytick", []);      # kill clrbar ticks too
    # (Some versions octave broken and need to kill ticks too to kill labels.)
    print("traveltime-$1-xyz-g2.0-nr.png", "-r200");  # PNG version
    #set(hRID, "color", [1 1 1]);  # Invisible but selectable in pdf version
    #set(findobj("tag", "OLCapt3"), "color", [1 1 1]); # Same for NR label
    newstr=get(findobj("tag", "OLCapt2"), "string");
    newstr=sprintf("%s\n%s","Channel:  Ex+Ey+Ez",newstr); # Add to OL2 label
    set(findobj("tag", "OLCapt2"), "string", newstr);
    #title("");  # Kill title
    print("tmp-traveltime-$1-xyz-g2.0-nr.pdf");       # PDF version
EOF
    pdfcroptobb "tmp-traveltime-$1-xyz-g2.0-nr.pdf" \
                "traveltime-$1-xyz-g2.0-nr.pdf"  # Crop to bb
}
# TT Curves:
produce_ttcurves WUS   0 159 8 0    # Station WUS
produce_ttcurves MAK 160 319 8 0    # Station MAK

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
