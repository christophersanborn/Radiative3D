#!/bin/bash
#
# Creates a subdirectory called wiki/ and populates it with stuff that
# I usually like to upload to the wiki.  Images are softlinked rather
# than copied as a space-saving measure and a hedge against
# re-plotting.  Filenames follow a defined convention.
#
# Also combines log files and such into wikitext suitable for the info
# pages that accompany images in the wiki.
#
# Run this script from inside a data directory created by a do-
# script.
#
#

if [ ! -r seis_000.octv ] && [ ! -r seisfiles/seis_000.octv ] ; then
    echo "Does not appear to be a Radiative3D output directory. Exiting."
    exit
fi

WIKIDIR=wiki

if [ -r wikify.incl ] ; then
    source wikify.incl   # Include file may have been written by do-script.
fi

# Create Wiki Dir:
if [ ! -d "$WIKIDIR" ]; then
    mkdir "$WIKIDIR"
fi

# Write a little renamer script:
cat > "$WIKIDIR/rename-files.sh" <<"EOF"
#!/bin/bash
if [ "$#" -ge 1 ]; then
  for file in *.{png,txt}; do 
    mv "$file" "$1 $file"
  done
fi
EOF
chmod u+x "$WIKIDIR/rename-files.sh"

# Figure out name of our do-script:
SHFILE=/dev/null
for file in *.sh; do
    if grep -qx '##  "DO"-script for a Radiative3D run:' "$file"
    then                # Catch ONLY file with R3D do marker
        SHFILE="$file"  # Also catches no more than one do- file
    fi                  #
done
TAG=${SHFILE%.sh}
TAG=${TAG#do-}

#
# "Copy" (actually link) the images that we typically keep for the wiki:
#
[ -z "$STA" ] && STA="WUS"    # If not already defined in wikify.incl
[ -z "$STB" ] && STB="MAK"    # then use these default values
[ -z "$STC" ] && STC="STC"    #
[ -z "$MODELCODE" ] && MODELCODE="X"  # E.g. "NSCP01" 

function LambdaLinkIf {  # Used on next few lines to link files into
                         # WIKIDIR but only if they exist.
  [ -r "$1" ] && ln -sf "../$1" "$WIKIDIR/$TAG $2"
}

LambdaLinkIf "seis_079.png" "Env079 $STA.png"
LambdaLinkIf "seis_119.png" "Env119 $STA.png"
LambdaLinkIf "seis_159.png" "Env $STA.png"
LambdaLinkIf "seis_239.png" "Env239 $STB.png"
LambdaLinkIf "seis_279.png" "Env279 $STB.png"
LambdaLinkIf "seis_319.png" "Env $STB.png"
LambdaLinkIf "seis_399.png" "Env399 $STC.png"
LambdaLinkIf "seis_439.png" "Env439 $STC.png"
LambdaLinkIf "seis_479.png" "Env $STC.png"
LambdaLinkIf "profile-LOP-WUS.png" "Model LOP-WUS.png"
LambdaLinkIf "profile-LOP-MAK.png" "Model LOP-MAK.png"
LambdaLinkIf "profile-WUS-MAK.png" "Model WUS-MAK.png"
LambdaLinkIf "profile-model.png" "Model $MODELCODE Elev.png"
LambdaLinkIf "seismap.png" "Seismometer Map.png"
for file in traveltime-???-xyz-*.png ; do
    ln -sf "../$file" "$WIKIDIR/$TAG TT ${file#traveltime-}"
done


#
# Now build the wikitext file:
#

LOGFILE=logfile
STDOUTFILE=stdout.txt
SVNFILE=svn-info.txt
MPARAMS=out_mparams.octv
OFILE="$WIKIDIR/$TAG Wikitext.txt"

if [ -r "$MPARAMS" ]; then
    numphonons=`grep NumPhonons -A 2 out_mparams.octv | tail -n 1`
fi

> "$OFILE"  # Clear file
echo "From data directory: " $(basename "`pwd`") >> "$OFILE"
echo  >> "$OFILE"
if [ -r description.txt ]; then
    cat description.txt >> "$OFILE"
    echo  >> "$OFILE"
fi
echo "Number of event phonons cast: " $numphonons >> "$OFILE"
echo  >> "$OFILE"
echo "Name of script file: [[#Shell File|$SHFILE]]" >> "$OFILE"
echo  >> "$OFILE"

echo "==== Log File ====" >> "$OFILE"
echo "<pre>" >> "$OFILE"
cat $LOGFILE >> "$OFILE"
echo "</pre>" >> "$OFILE"
echo  >> "$OFILE"

echo "==== Output ====" >> "$OFILE"
echo "<pre>" >> "$OFILE"
cat $STDOUTFILE >> "$OFILE"
echo "</pre>" >> "$OFILE"
echo  >> "$OFILE"

echo "==== Shell File ====" >> "$OFILE"
echo '<syntaxhighlight lang="bash">' >> "$OFILE"
cat $SHFILE >> "$OFILE"
echo "</syntaxhighlight>" >> "$OFILE"
echo  >> "$OFILE"

if [ -r do-figsonly.sh ]; then
  echo "==== Figure File ====" >> "$OFILE"
  echo '<syntaxhighlight lang="bash">' >> "$OFILE"
  cat do-figsonly.sh >> "$OFILE"
  echo "</syntaxhighlight>" >> "$OFILE"
  echo  >> "$OFILE"
fi

if [ -r do-normcompare.sh ]; then
  echo "==== Normcurve File ====" >> "$OFILE"
  echo '<syntaxhighlight lang="bash">' >> "$OFILE"
  cat do-normcompare.sh >> "$OFILE"
  echo "</syntaxhighlight>" >> "$OFILE"
  echo  >> "$OFILE"
fi

echo "==== Svn File ====" >> "$OFILE"
echo "<pre>" >> "$OFILE"
cat $SVNFILE >> "$OFILE"
echo "</pre>" >> "$OFILE"
echo  >> "$OFILE"


#
# Tar up seis_files:
#
#tar -czf seis_nnn.tar.gz seis_[0-9][0-9][0-9].octv

#
# Now make a nice little tarball for easy downloading:
#

TARFILE="$WIKIDIR"/"WikiDir-$TAG.tar.gz"
touch "$TARFILE"  # Prevents warning about $WIKIDIR modified as tar is
                  # running.

tar -czhf "$TARFILE" \
    --transform="s,$WIKIDIR,$TAG," \
    --exclude='*.tar.gz' \
    "$WIKIDIR"/
