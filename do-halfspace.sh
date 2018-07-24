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
##  This DO-SCRIPT sets up a run based on a HALFSPACE model with up to
##  ONE interface layer.
##

## One-liner description: (Keep this BRIEF.)
##
INTENT="Example halfspace simulation."
CAMPAIGN="" ## (e.g. "For AGU Poster Dec 2016" or some such.)

SIMTARGET="waveform"          # Choice: 'waveform' or 'video'. Affects
                              # defaults not otherwise specified.
MODIDX=40  # HALFSPACE        # Model Index: Selects from custom coded models.
                              # 1: Lop Nor (base or moho depends on COMPARGS),
                              # 5: North Sea Crust Pinch model,
                              # 8: Crust Upthrust model
                              # 40: Halfspace model

event=eq        # 'eq' or 'expl' - See case statement below for choices
EQISOFRAC=0.0   # Iso fraction for EQ event (choose from range [-1.0, 1.0])
                # (Or choose as Iso Angle in range [-90, -1), (1, 90].)

##
##  Epicenter and Station Locations
##

EPIXYZ=0,0,0                  #
STAXYZ=183.85,183.85,0        # Station 'A', 260 km to the NE; Pg azimuth
STBXYZ=260,0,0                # Station 'B', 260 km to the East; Lg azimuth
STCXYZ=240.21,-99.5,0         # Station 'C', 260 km to the ESE; mixed Pg/Lg

FREQ=2.0                      # Phonon frequency to model
NUMPHONS=10M                  # Number of phonons to spray (Recommend: 140M)
RECTIME=200                   # Recording duration of seismometers (seconds).
BINSIZE=0.50                  # Seismometer time-bin size in seconds
GATHER=0.105,10.0             # Initial,Terminal gather radius, in kilometers.
R0=2.737                      # Range of first receiver in array
NSEIS=48                      # Number of receivers per linear array
CYLRAD=900                    # Total (cylindrical) radius of model.
                              # (Phonons exceeding this radius from XY
                              # = (0,0), or this travel time from t=0,
                              # are abandoned.)

SCAT1=0.8,0.01,1.0,0.5,1000   # Scat Args (nu,eps,a,kappa) and Q in top layer
SCAT2=0.8,0.01,1.0,0.5,1000   # Scat Args (nu,eps,a,kappa) and Q in bot layer
VPVS1=6.40,3.63,2.83,-60      # Seismic attrs (Vp, Vs, rho, z-bottom) top layer
VPVS2=6.40,3.63,2.83,-400     # Seismic attrs (Vp, Vs, rho, z-bottom) bot layer
                              # (Note: Q's specified are Q_s values.
                              # Q_p is computed from assumption that
                              # Q_kappa is infinite.)
COMPARGS=$SCAT1,$SCAT2,$VPVS1,$VPVS2

#FLATTEN="--flatten"          # If set to "--flatten", apply Earth-flattening
FLATTEN=""                    # transformation to depths and velocities.
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
        SOURCELOC=0,0,-5      # Origin, 5km deep
        SOURCETYP=EXPL        #
        ;;
    eq)     # Generic earthquake
        SOURCELOC=0,0,-5      # Origin, 5km deep
        SOURCETYP=SDR,0,90,0,$EQISOFRAC     # Basic Strike-Slip
        ;;
    *)
        echo Event code not recognized.
        exit
        ;;
esac

SEISORIG1=$EPIXYZ             # Array locations:
SEISORIG2=$EPIXYZ             #
SEISORIG3=$EPIXYZ             #
SEISDEST1=$STAXYZ             #
SEISDEST2=$STBXYZ             #
SEISDEST3=$STCXYZ             #
SEIS1=--seis-p2p=$SEISORIG1,$SEISDEST1,$R0,$GATHER,$NSEIS
SEIS2=--seis-p2p=$SEISORIG2,$SEISDEST2,$R0,$GATHER,$NSEIS
SEIS3=--seis-p2p=$SEISORIG3,$SEISDEST3,$R0,$GATHER,$NSEIS


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
  STA="PGX"
  STB="LGX"
  STC="PLG"
  MODELCODE="HLF"
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
  ## Seismometer Map:
  seismap("seisfiles/seis_???.octv","out_mparams.octv",
          [-50 274 -47 196], [45 95 118], {"PGX","LGX","PLG"});
  legend("hide"); # (We're not using separate P and S gather anymore)
  print("seismap.png","-r200");   # PNG version for wiki
  print("tmp-seismap.pdf");       # PDF version for publication
  ## LOP to MAK:
  modelplot("griddump.txt", [1 2], {"SRC","PGX"},
            "out_mparams.octv",
            "seisfiles/seis_%03d.octv",
            [0,1:2:47]);
  print("profile-LOP-MAK.png","-r200"); # PNG version for wiki
  #delete(findall("tag","sourcemark")); # Remove source marker for publication
  #delete(findall("tag","seismometer"));# Remove seismometer marks if desired
  print("tmp-profile-LOP-MAK.pdf");     # PDF version for publication
  ## LOP to WUS:
  modelplot("griddump.txt", [1 3], {"SRC","LGX"},
            "out_mparams.octv",
            "seisfiles/seis_%03d.octv",
            [48,49:2:95]);
  print("profile-LOP-WUS.png","-r200");
  ## WUS to MAK:
  modelplot("griddump.txt", [3 2], {"LGX","PGX"},
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
for file in seisfiles/{seis_{000,048,096},seis_?[02479][357]}.octv ; do
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
for file in seisfiles/{seis_047.octv,seis_095.octv,seis_143.octv} ; do
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
    NORMCURVE = arrayimage(AR, [1 1 1], 2.0, [0 200]);
    save("traveltime-$1-xyz-normcurve.octv", "NORMCURVE");
    NC2 = NORMCURVE;
    NORMCURVE = normcurve_fitpowerlaw(NORMCURVE,4,48);
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
    arrayimage(AR, [1 1 1], 2.0, [0 200], NORMCURVE, $4);
    annotate_array([6.40,3.63],{"Pg","Lg"},230); # (mark phase velocities)
    hRID = papertext(0.03, 0.00, "$RUNID",   ## Plot RunID in lowerleft corner
              "verticalalignment", "bottom", "horizontalalignment", "left",
              "fontsize", 5, "fontweight", "demi");
    print("traveltime-$1-xyz-g2.0.png", "-r200");
    ## VARIANT 2: And make one where norm relative to local peak and average:
    arrayimage(AR, [1 1 1], 2.0, [0 200], 0.5, ${5:-0});
    annotate_array([6.40,3.63],{"Pg","Lg"},230); # (mark phase velocities)
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
produce_ttcurves PGX   0  47 8 0    # Station LGX
produce_ttcurves LGX  48  95 8 0    # Station PGX
produce_ttcurves PLG  96 143 8 0    # Station PLG

# Lapse Time Curves (Template Function):
produce_lapsecurves() {  # $1:     filename tag (eg station code if relevant)
                         # $2, $3: seis indices
                         # $4, $5: lapse times (default 0, 30)
  phaseedge="[3.63, 0]"
  axesxyz="[1,1,1]"
  axestxt="xyz"
  geospread=2
  LT1="${4:-0}"
  LT2="${5:-30}"
  QIntS="$(grep -A 10 '#  R3D_GRID:' stdout.txt | grep '^  0   0   0' | awk '{ print $11 }')"
  QScatS="$(grep -A 6 '#  BEGIN SCATTERER DUMP:' stdout.txt | grep -m 1 '^ ' | awk '{ print $5 * $9 }')"
  echo "$QScatS $QIntS" > QScatIntS.txt
  octave -qf <<EOF

    AR1 = array("seisfiles/seis_%03d.octv",$2,$3);
    ARs = {AR1}; nAR = 1;
    QScatInt = {load("QScatIntS.txt")};
    otherdirs = {"other"};
    for (i=1:25); otherdirs{end+1} = sprintf("other%d",i); end
    for (i=1:length(otherdirs))
      if (exist(otherdirs{i},"dir")==7)
        filepattern = sprintf("%s%s",otherdirs{i},"/seisfiles/seis_%03d.octv");
        AR = array(filepattern,$2,$3);
        ARs{end+1} = AR; nAR += 1;
        QScatInt{end+1} = load(sprintf("%s/QScatIntS.txt",otherdirs{i}));
      end
    end
    LineWidths = {7,6,5,4,3,2,1};
    LineStyles = {"-", "-.", ":"};

    LgEdge = $phaseedge;  # [v, t0]
    LWindow1 = $LT1 + [0 15];
    LWindow2 = $LT2 + [0 70];

    lapsetimebaseplot();
    global linestyle;

    R1R2data = [];
    for i=1:nAR
      linestyle.style = LineStyles{mod(i-1,length(LineStyles))+1};
      linestyle.width = LineWidths{mod(i-1,length(LineWidths))+1};
      data = lapsetimecurve(ARs{i}, QScatInt{i}, $axesxyz,
                            $geospread, LgEdge, LWindow1, LWindow2);
      R1R2data = [R1R2data; data];
    end
    R1R2data

    axis([0 250, 0.0001 10]);
    titletxt = sprintf(
        "Blue: Early window (%d to %d s); Green: Late window (%d to %d s)",
        LWindow1(1), LWindow1(2), LWindow2(1), LWindow2(2)
    );
    title(titletxt, "fontweight", "bold");

    hlgnd = findobj(gcf(),"type","axes","Tag","legend");
    set(hlgnd, "fontweight", "bold", "FontSize", 8);

    print("tmp-LapseTime-$1-$axestxt.pdf");
    print("LapseTime-$1-$axestxt.png","-r200");

    # Now plot the R1 R2 data:
    fehlerR1R2plot(R1R2data);
    print("tmp-FehlerR1R2-$1-$axestxt.pdf");

EOF
  pdfcroptobb "tmp-LapseTime-$1-$axestxt.pdf" \
              "LapseTime-$1-$axestxt.pdf"  # Crop to bb
  pdfcroptobb "tmp-FehlerR1R2-$1-$axestxt.pdf" \
              "FehlerR1R2-$1-$axestxt.pdf"  # Crop to bb
}
# Lapse Time curves:
produce_lapsecurves LGX 48 95 0 30      # Lg Lapse Time curve at LGX

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
