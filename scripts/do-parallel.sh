#!/bin/bash
##
## Start a parallel batch of do-sequential scripts each in their own
## individual screen session.
##
## If optional arguments $1 and $2 are provided, they are interpreted
## respectively as the minimum and maximum index number of the do-scripts to
## run.  This can be used to selectively re-run certain ranges of scripts
## after the initial run, or to run the batch set in separate chunks,
## etc. (The arguments are passed to the do-sequential files and are actually
## handled there.)
##
## The do-sequential scripts are used with large batches of do-scripts to
## apportion them across available CPU cores.  The do-scripts and
## do-sequential scripts are pregenerated by a dot-script (do-template
## script) to genrate a large set of simulation conditions and determine how
## many sequential streams to apportion across.
##
TAG=`date +"%Y%m%d-%H%M%S" | shasum`
TAG="P-${TAG:0:6}"
count=0

for file in do-sequential-[0-9]*.sh; do
    idx=${file%.sh}
    idx=${idx##*-}
    scrsesh="$TAG-$idx-R3D"
    screen -S $scrsesh -dm sh -c "./$file $*"
    ((count+=1))
done
echo "Started $count screen sessions running with prefix $TAG."
