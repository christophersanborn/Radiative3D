#!/bin/bash
# sprayhisto.sh
# Plots a histogram of the spray angles from from GEN events to give a
# visualization of the directional dependnece of P and S radiation.
# Intended to test manipulation of moment tensors.
#
RAD3D=../../main

REPORTS=reports.txt
GENPDATA=gendataP.txt  # Two columns: Theta and Phi (w.r.t. ENU basis)
GENSDATA=gendataS.txt  #
NUMPH=36000

SUFFIX="Before"

function sprayhisto_main {

  EVENT=SDR,0,90,0
  runandplot $EVENT "$EVENT ($SUFFIX Fix)" "_$SUFFIX"

  EVENT=SDR,22.5,90,0
  runandplot $EVENT "$EVENT ($SUFFIX Fix)" "_$SUFFIX"

  EVENT=SDR,45,90,0
  runandplot $EVENT "$EVENT ($SUFFIX Fix)" "_$SUFFIX"

  EVENT=USGS,1,0,0,0,0,0      ## RR,tt,pp,...
  runandplot $EVENT "$EVENT (LVD UD) ($SUFFIX)" "_$SUFFIX"

  EVENT=USGS,0,1,0,0,0,0      ## rr,TT,pp,...
  runandplot $EVENT "$EVENT (LVD NS) ($SUFFIX)" "_$SUFFIX"

  EVENT=USGS,0,0,1,0,0,0      ## rr,tt,PP,...
  runandplot $EVENT "$EVENT (LVD EW) ($SUFFIX)" "_$SUFFIX"

}

function runandplot {
  event="$1"
  title="$2"
  suffix="$3"

  $RAD3D --reports=ALL_ON --source=$event -N $NUMPH --report-file=$REPORTS

  cat $REPORTS | grep "^GEN:" | grep "  P  " | \
      awk '{ print $14, $15 }' > $GENPDATA
  cat $REPORTS | grep "^GEN:" | grep "  S  " | \
      awk '{ print $14, $15 }' > $GENSDATA
  rm $REPORTS

  octave -qf <<EOF
    ThetaPhiENUP = load("$GENPDATA");
    ThetaPhiENUS = load("$GENSDATA");
    sprayhisto(ThetaPhiENUP,ThetaPhiENUS,"$title");
    print("SprayHist_${event}${suffix}.png", "-r200");
EOF
}

sprayhisto_main
