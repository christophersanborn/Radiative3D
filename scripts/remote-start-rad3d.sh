#!/bin/bash
#
#  Starts one or more do-script(s) in screen sessions on a remote machine
#
#  Takes same syntax as remote-inst-rad3D.
# 
#  Depends on remote-inst-rad3d.sh to push script files to remote machine
#
DEST=/scratch/`whoami`/WorkingCopies/Radiative3D
RINSTALLER=`dirname $0`/remote-inst-rad3d.sh

# Expects two or more arguments
if [ $# -lt 2 ]; then
    echo "Usage:  $0 machine do-script.sh ..."
    echo "Starts do-script(s) on remote machine. Ensures R3D is"
    echo "installed and up-to-date first."
    exit
fi

HOST=$1   # Get remote host from first arg
shift     # Kill first arg and shift remainder down

# Install or update and push files first:
$RINSTALLER $HOST $@

# Remote start scripts in nice little screensessions:
for file in "$@"; do 
    echo "Starting `basename $file`"
    tag=`basename $file`
    tag=${tag%.*}
    tag=${tag#do-}
    ssh -T $HOST /bin/bash <<____EOF
        screen -S `basename ${file%.*}` -dm \
        sh -c "cd $DEST; ./`basename $file` $tag; echo press any key; read"
____EOF
    sleep 2   # Guarantees different timestamps
done
