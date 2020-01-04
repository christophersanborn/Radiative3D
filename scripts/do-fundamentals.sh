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
## Helpers:
##
function SecondsToHours {
    printf %.2f\\n "$((10**9 * $1/3600))e-9"
}
function SecondsDeltaToHours {
    local Delta=$(($1 - $2))
    echo $(SecondsToHours $Delta)
}
function VCSVersion {
    if [ -d .svn ]; then
        svnversion
    elif [ -d .git ]; then
        git rev-parse --short HEAD
    else
        echo "Unversioned directory"
    fi
}
function VCSStatusInfo {
    if [ -d .svn ]; then
        echo "Revision number (svnversion): $(svnversion)"
        svn info
        svn st
        svn diff
    elif [ -d .git ]; then
        echo "Revision number (git): $(git rev-parse --short HEAD)"
        git status
        git diff
    else
        echo "Unversioned directory"
    fi
}
function RealPath {
    # realpath not a base install on Mac (although can be added with
    # `brew install coreutils`). We fall back to just echoing the path if
    # not installed.
    if which realpath > /dev/null; then
        realpath "$1"
    else
        echo "$1"
    fi
}

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
        template)
            [ -z "$FINTAG" ]     && FINTAG=R3BATCH
            # No further defualts if we're just templating a batch.
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
## FUNCTION:  CreateBatchDirectory()
##
##   Sets up a top-level output directory for batch output.  Directory
##   will be named with a timestamp and an extra token if $1 is set.
##   Useage:
##
##     CreateBatchDirectory $@
##
##   (The $@ is needed to send script args to this function, the main
##   one is a tag id in the output directory name.)
##
##   @sa CreateBatchSubDirectory().
##   @sa CreateOutputDirectory().
##
##   ASSUMPTIONS:
##
##   Function assumes that $OUT_BASE and $R3D_EXE are already defined.
##
##   SIDE EFFECTS:
##
##   A directory will be created and populated. $R3D_EXE will be modified.
##   $outdir will be set, along with $revdir and $doshdir for revision state
##   and do-script library.  $TMPL_SRC will point to the extracted template
##   source file.  $RUNID will be set.  A logfile will be started and
##   $LOGFILE will be set.
##
function CreateBatchDirectory {

    local SHFILE="$(basename "$0")"     # Guess a default RunID from
    RUNID=${SHFILE%.sh}                 # the shell file name,
    RUNID=${RUNID#dot-}                 #
    [ -z "$1" ] || RUNID="$1"           # but use $1 if it was provided.

    [ -z "$FINTAG" ] && FINTAG=R3DBATCH
    local tstamp=`date +"%Y%m%d-%H%M%S"`
    local outdirname="$tstamp"-"$RUNID"-"$FINTAG"
    outdir="$OUT_BASE/$outdirname"
    revdir="$outdir/revision-files"         # Retain revision used for sim batch
    doshdir="$outdir/do-scripts"            # Put do-scripts in here
    combodir="$outdir/cumulative"
    sumfigdir="$outdir/summaryfigs"
    mkdir -p "$sumfigdir"
    mkdir -p "$combodir"
    mkdir -p "$revdir"
    mkdir -p "$doshdir"
    echo "Summation of repeated runs goes in here." > "$combodir/WhatsThis.txt"
    echo "Do-scripts in here but RUN them from dir above." > "$doshdir/WhatsThis.txt"

    if [ -L "$OUT_BASE/previous" ]; then    # Generate "latest" and "previous"
        rm "$OUT_BASE/previous"             # links to make it easier to cd
    fi                                      # into new output
    if [ -L "$OUT_BASE/latest" ]; then
        mv "$OUT_BASE/latest" "$OUT_BASE/previous"
    fi
    ln -sf "$outdirname" "$OUT_BASE/latest"
    #

    cp "$0" "$revdir"/                # Copy do-script
    cp $R3D_EXE "$revdir"/            # Keep the exe for replicability
    R3D_EXE="$revdir"/$(basename $R3D_EXE)  # Run the copy instead of the orig
    echo "(BRIEF description of batch goes here.)" > "$outdir"/description.txt
    [ -z "$INTENT" ] || echo "$INTENT" >> "$outdir"/description.txt
    [ -z "$CAMPAIGN" ] || echo "$CAMPAIGN" >> "$outdir"/description.txt
    VCSStatusInfo > "$outdir"/svn-info.txt # Record detailed revision info
    cp vis/seisplot/*.m "$revdir"/          # Used in post-process
    cp vis/scattervid/*.m "$revdir"/        #      ''
    cp vis/scattervid/preprocess.sh "$revdir"/ #   ''
    cp scripts/do-fundamentals.sh "$revdir"/   #   ''
    cp scripts/combineresults.sh "$revdir"/ #      ''
    cp scripts/wikify.sh "$revdir"/         #      ''
    ln -sr "$combodir" "$sumfigdir/datasrc"
    ln -sr "$revdir"/*.m "$sumfigdir"/

    LOGFILE="$outdir"/logfile
    > "$LOGFILE"
    echo "Machine: " `hostname` >> "$LOGFILE"
    echo "User: " `whoami` >> "$LOGFILE"
    echo "Revision: " `VCSVersion` >> "$LOGFILE"

    # Carve-off and save the template part so we can stream-edit it after:
    TMPL_SRC="$outdir"/"do-template-src"
    cat "$0" | grep -A 10000 -B 1 -E "^## ___TEMPLATE_START___$" | \
               grep -A 0 -B 10000 -E "^## ___TEMPLATE_END___$" | \
               sed 's/ ___TEMPLATE_START___//' \
               > "$TMPL_SRC"
    cat "$0" | grep -A 10000 -B 1 -E "^## ___SUMMARY_FIGS_START___$" | \
               grep -A 0 -B 10000 -E "^## ___SUMMARY_FIGS_END___$" | \
               sed 's/ ___SUMMARY_FIGS_START___//' \
                   > "$sumfigdir/do-summaryfigs.sh"
    chmod a+x "$sumfigdir/do-summaryfigs.sh"
}

function CreateBatchSubDirectory {

    local SHFILE="$(basename "$0")"     # Guess a default RunID from
    local RUNID=${SHFILE%.sh}           # the shell file name,
          RUNID=${RUNID#do-}            #
    [ -z "$1" ] || RUNID="$1"           # but use $1 if it was provided.

    [ -z "$FINTAG" ] && FINTAG="R3D"
    local tstamp="0"                    # TODO: increment if already used
    local outdirname="runs/$tstamp"-"$RUNID"-"$FINTAG"
    while [ -e "$outdirname" ]; do
        ((tstamp++))
        outdirname="runs/$tstamp"-"$RUNID"-"$FINTAG"
    done
    outdir="$outdirname"                # No longer has a base dir; wrt current dir
    revdir="revision-files"             # Where revision files are stored
    combodir="cumulative"
    mkdir -p "$outdir"
    touch "$outdir.inprogress"

    cp "$0" "$outdir"/                  # Copy do-script (Used only to get RUNID)
    R3D_EXE="$revdir"/$(basename $R3D_EXE)  # Run the revision copy
    echo "(BRIEF description of batch goes here.)" > "$outdir"/description.txt
    [ -z "$INTENT" ] || echo "$INTENT" >> "$outdir"/description.txt
    [ -z "$CAMPAIGN" ] || echo "$CAMPAIGN" >> "$outdir"/description.txt
    ln -sr "$revdir"/*.m "$outdir"/             # Used in post-process
    ln -sr "$revdir"/preprocess.sh "$outdir"/   #   ''

    LOGFILE="$outdir"/logfile
    > "$LOGFILE"
    echo "Machine: " `hostname` >> "$LOGFILE"
    echo "User: " `whoami` >> "$LOGFILE"
    echo "Revision: " `VCSVersion` >> "$LOGFILE"

    # Carve-off and save figure-generation part of do-script in case
    # we want to re-run the figures later (as might happen if we
    # improve or tweak the plotting routines after a run).
    cat "$0" | grep -A 10000 -B 1 -E "^## ___FIG_GEN_START___$" | \
               grep -A 0 -B 10000 -E "^## ___FIG_GEN_END___$" \
               > "$outdir"/do-figsonly.sh
    chmod u+x "$outdir"/do-figsonly.sh

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
    echo `hostname`:`RealPath "$outdir"` >> $rcptfile

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
    VCSStatusInfo > "$outdir"/svn-info.txt # Record detailed revision info
    cp vis/seisplot/*.m "$outdir"/          # Used in post-process
    cp vis/scattervid/*.m "$outdir"/        #      ''
    cp vis/scattervid/preprocess.sh "$outdir"/ #   ''
    cp scripts/wikify.sh "$outdir"/         #      ''

    LOGFILE="$outdir"/logfile
    > "$LOGFILE"
    echo "Machine: " `hostname` >> "$LOGFILE"
    echo "User: " `whoami` >> "$LOGFILE"
    echo "Revision: " `VCSVersion` >> "$LOGFILE"

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
    local TIMEELAPSED=$(SecondsDeltaToHours $TIMEEND $TIMEBEGIN)
    echo "($TIMEELAPSED Hours Run-Time)" >> "$LOGFILE"

}

######
## FUNCTION:  SpiralIndex n=$1 i0=$2 j0=$3
##
##  Given an index n as a counter in an i,j index space, return "i j" after
##  tracing a spiral path from centerpoint index pair (i0,j0). This is used
##  to iterate a 2-D parameter space "from the center outward" so that we
##  can view early results clustered around the center and the later results
##  will paint the perimeter.  n=0 returns "i0 j0". n>0 traces path like:
##
##      13 - 14 - 15 ...         +i -->
##       |                    +j
##      12    2 -- 3 -- 4      |
##       |    |         |     \ /
##      11    1 -- 0    5      V
##       |              |
##      10 -- 8 -- 7 -- 6
##
function SpiralIndex {

    local n="$1"
    [ "$n" -ge "0" ] || return  # Integer check

    local i="$2"
    local j="$3"
    [ "$i" -ge "0" ] || return  # Integer check
    [ "$j" -ge "0" ] || return  # Integer check

    if [ "$n" -eq "0" ]; then
        echo "$i $j"
        return
    fi

    local p=1
    local step=-1
    local idx=0
    while [ "$p" -lt "40" ]; do  # (Set large sanity limit)
        local q=1
        local qmax="$p"; ((qmax*=2))
        while [ "$q" -le "$qmax" ]; do
            if [ "$q" -le "$p" ]; then
                ((j+=step))
            else
                ((i+=step))
            fi
            ((idx+=1))
            if [ "$idx" -eq "$n" ]; then
                echo $i $j
                return
            fi
            ((q+=1))
        done
        ((p+=1))
        ((step*=-1))
    done
    # shouldn't get here unless n too big

}
