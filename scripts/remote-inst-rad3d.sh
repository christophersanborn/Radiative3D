#!/bin/bash
#
#   Install Radiative3D and copy do-script files to a remote host.
#   E.g.,
#
#     ./remote-inst-rad3d.sh loki do-loki-run1.sh do-loki-run2.sh
#
#   will put Radiative3D and two do-scripts into the destination
#   folder on loki.
#
#
DEST=/scratch/`whoami`/WorkingCopies/Radiative3D
REPO=https://rainbow.phys.uconn.edu/svn/Radiative3D/trunk

if [ $# -lt 1 ]; then
    echo "Usage:   $0  [user@]host  [file]  ..."
    echo "Installs Radiative3D on remote host and copies file(s) to the install directory."
    exit
fi

HOST=$1   # Get remote host from first arg
shift     # Kill first arg and shift remainder down

# Install Radiative3D:
echo "INSTALLING Radiative3D on $HOST; You may need to type your password:"
ssh -T $HOST /bin/bash <<EOF
  mkdir -p $DEST
  cd $DEST
  if [ ! -f $DEST/main.cpp ]; then 
    svn co --non-interactive --trust-server-cert $REPO ./
  else
    svn up --non-interactive --trust-server-cert
  fi
  make
  pwd
EOF

# Now copy files to dest directory:
if [ $# -ge 1 ]; then
    echo "COPYING files $@ to $HOST; You may need to type your password:"
    scp "$@" $HOST:$DEST/
fi
