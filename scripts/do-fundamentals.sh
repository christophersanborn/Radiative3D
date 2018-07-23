#!/bin/bash
##
##  DO-SCRIPT FUNDAMENTALS
##
##  (Library of routines for Radiative3D Do-Scripts.)
##
##  Include (source) this file in your own do-scripts, with a line like this:
##
##    source scripts/do-fundamentals.sh
##
##

 
######
## FUNCTION:  PopDefaults()
##
##   Populates default values for control veriables that are not
##   already specified, based on a selected "target" simulation type.
##   E.g., certain command line option and their values are determined
##   mainly by whether you are producing envelopes or videos, but do
##   not otherwise change from one simulation to another. These are
##   packed into here and chosen via the $1 argument.
##
##   SIDE AFFECTS:
##
##   Sets (but will not override if already exists) $R3D_EXE,
##   $OUT_BASE, $TOA_DEGREE, $REPORTS, (and others, tbd)
##
##
function PopDefaults {

    [ -z "$R3D_EXE" ]  && R3D_EXE=./main  # Executable to use
    [ -z "$OUT_BASE" ] && OUT_BASE=data   # Parent directory in which
                                          # to put output
                                          # subdirectories

    case "$1" in
        waveform)
            [ -z "$TOA_DEGREE" ] && TOA_DEGREE=9
            [ -z "$REPORTS" ]    && REPORTS=ALL_OFF
            [ -z "$FINTAG" ]     && FINTAG=R3D
            ;;
        video)
            [ -z "$TOA_DEGREE" ] && TOA_DEGREE=8
            [ -z "$REPORTS" ]    && REPORTS=ALL_ON
            [ -z "$FINTAG" ]     && FINTAG=R3VID
            ;;
        *)
            echo Simulation target not recognized.
            exit
            ;;
    esac

}


######
## FUNCTION:  EventParams()
##
##   A selector for named event types.  First arg is event tag.
##   Subsequent args depend on selection.
##
##   USAGE:
##
##     SOURCETYP=$(EventParams [args])
##
##
function EventParams {

    local result="EQ"
    local isofrac=0

    case "$1" in
        expl)   # Generic explosion
            SOURCETYP=EXPL                    #
            ;;
        Selby|selby)
            if [ $# -gt 1 ]; then
                isofrac=$2
            fi
            result=SDR,125,40,90,$isofrac
            ;;
        *)
            echo Event code not recognized.
            exit
            ;;
    esac
    echo "$result"
}


######
## FUNCTION:  CheckBuildStatus()
##
##   Checks whether Make indicates that $R3D_EXE is up-to-date.  Exits
##   script with a message if it is not.
##
function CheckBuildStatus {
    if ! make -q
    then
        echo Make indicates $R3D_EXE is not up to date.
        echo Please run make and re-run this script.
        exit
    fi
}


######
## FUNCTION:  CreateOutputDirectory()
##
##   Sets up directory in which to put simulation output, graphical
##   results, and and in which to store input records and scripts,
##   etc.  Directory will be named with a timestamp and an extra token
##   if $1 is set.  Useage:
##
##     CreateOutputDirectory $@
##
##   (The $@ is needed to send script args to this function, the main
##   one is a tag id in the output directory name.)
##
##   ASSUMPTIONS:
##
##   Function assumes that $OUT_BASE and $R3D_EXE are already defined.
##
##   SIDE EFFECTS:
##
##   A directory will be created and populated. $R3D_EXE will be
##   modified.  $outdir will be set.  A logfile will be started and
##   $LOGFILE will be set.
##
function CreateOutputDirectory {

    local SHFILE="$(basename "$0")"     # Guess a default RunID from
    local RUNID=${SHFILE%.sh}           # the shell file name,
          RUNID=${RUNID#do-}            #
    [ -z "$1" ] || RUNID="$1"           # but use $1 if it was provided.

    [ -z "$FINTAG" ] && FINTAG=R3D
    local tstamp=`date +"%Y%m%d-%H%M%S"`
    local outdirname="$tstamp"-"$RUNID"-"$FINTAG"
    outdir="$OUT_BASE/$outdirname"
    mkdir -p "$outdir"
    local rcptfile="`basename ${0%.sh}`.$$.runrcpt"
    echo `hostname`:`realpath "$outdir"` >> $rcptfile

    if [ -L "$OUT_BASE/previous" ]; then    # Generate "latest" and "previous"
        rm "$OUT_BASE/previous"             # links to make it easier to cd
    fi                                      # into new output
    if [ -L "$OUT_BASE/latest" ]; then
        mv "$OUT_BASE/latest" "$OUT_BASE/previous"
    fi
    ln -sf "$outdirname" "$OUT_BASE/latest"
    #

    cp "$0" "$outdir"/                # Copy do-script
    cp $R3D_EXE "$outdir"/            # Keep the exe for replicability
    R3D_EXE="$outdir"/$(basename $R3D_EXE)  # Run the copy instead of the orig
    echo "(BRIEF description of run goes here.)" > "$outdir"/description.txt
    [ -z "$INTENT" ] || echo "$INTENT" >> "$outdir"/description.txt
    [ -z "$CAMPAIGN" ] || echo "$CAMPAIGN" >> "$outdir"/description.txt
    echo "Revision number (svnversion): " \
        `svnversion` > "$outdir"/svn-info.txt # Record detailed revision info
    svn info >>"$outdir"/svn-info.txt       #
    svn st >> "$outdir"/svn-info.txt        #
    svn diff >> "$outdir"/svn-info.txt 
    cp vis/seisplot/*.m "$outdir"/          # Used in post-process
    cp vis/scattervid/*.m "$outdir"/        #      ''
    cp vis/scattervid/preprocess.sh "$outdir"/ #   ''
    cp scripts/wikify.sh "$outdir"/         #      ''

    LOGFILE="$outdir"/logfile
    > "$LOGFILE"
    echo "Machine: " `hostname` >> "$LOGFILE"
    echo "User: " `whoami` >> "$LOGFILE"
    echo "Revision: " `svnversion` >> "$LOGFILE"

    # Carve-off and save figure-generation part of do-script in case
    # we want to re-run the figures later (as might happen if we
    # improve or tweak the plotting routines after a run).
    cat "$0" | grep -A 10000 -B 1 -E "^## ___FIG_GEN_START___$" | \
               grep -A 0 -B 10000 -E "^## ___FIG_GEN_END___$" \
               > "$outdir"/do-figsonly.sh
    chmod u+x "$outdir"/do-figsonly.sh

}


######
## FUNCTION:  RunSimulation()
##
##   Runs the actual Radiative3D simulation. Logs the begin time, run
##   time, and command line to the log file.
##
##   ASSUMPTIONS:
##
##   Assumes that $R3D_EXE and $LOGFILE have been defined and properly
##   set, and that all variables used in defining the command line
##   have been defined.  Typically, these are defined in the body of
##   the main do-script.
##
##   SIDE EFFECTS:
##
##   Begins simulation, appends to >>$LOGFILE
##
##
function RunSimulation {
    echo Begin Radiative3D Run: `date +"%Y.%m.%d-%H:%M:%S"` >> "$LOGFILE"
    local binsizearg=""
    [ -z "$BINSIZE" ] || binsizearg="--binsize=$BINSIZE"
    local cylradarg=""
    [ -z "$CYLRAD" ] || cylradarg="--range=$CYLRAD"
    local modelarg=""
    [ -z "$MODIDX" ] || modelarg="--grid-compiled=$MODIDX"
    local TIMEBEGIN=`date +"%s"`
    local R3D_CMDLN="\
$R3D_EXE --reports=$REPORTS \
        --output-dir=\"$outdir\" \
        --report-file=reports.dat \
        --mparams-outfile=out_mparams.octv \
        --num-phonons=$NUMPHONS \
        --toa-degree=$TOA_DEGREE \
        --source=$SOURCETYP \
        --source-loc=$SOURCELOC \
        --frequency=$FREQ \
        --timetolive=$RECTIME \
        $binsizearg \
        $modelarg \
        $cylradarg \
        $FLATTEN \
        $MFPOVERRIDE \
        $NODEFLECT \
        --model-args=$COMPARGS \
        --dump-grid \
        $SEIS1 \
        $SEIS2 \
        $SEIS3 \
        $ADDITIONAL \
"
    echo $R3D_CMDLN >> "$LOGFILE"
    eval "$R3D_CMDLN" | tee "$outdir"/stdout.txt
    echo End__ Radiative3D Run: `date +"%Y.%m.%d-%H:%M:%S"` >> "$LOGFILE"
    local TIMEEND=`date +"%s"`
    local TIMEELAPSED=`echo "scale=2; ($TIMEEND-$TIMEBEGIN)/3600" | bc`
    echo "($TIMEELAPSED Hours Run-Time)" >> "$LOGFILE"

}
