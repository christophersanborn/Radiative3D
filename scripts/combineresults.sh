#!/bin/bash
#
#  Combines results from two runs to yield a dataset with a higher
#  phonon count.  To combine three or more, need to run multiple
#  times.  E.g. 'combine A B AB; combine AB C ABC'
#
#  Note: Can't handle spaces in seis file name.  I think spaces in
#  directory names are OK if properly quoted.
#
seispattern='seisfiles/seis_???.octv'  # Breaks on spaces. Can't figure out how
                                       # to make spaces work with the globbing
                                       # below.
parampattern='out_mparams.octv' # Combine.m also works for parameter file,
                                # e.g. to compute summed phonon count.

# Take a couple guesses at location of combine.m
COMB_M="$(dirname "${BASH_SOURCE[0]}")"/../vis/seisplot/combine.m
[ -f "$COMB_M" ] || COMB_M="$(dirname "${BASH_SOURCE[0]}")"/combine.m

print_usage(){
    echo "  Usage:"
    echo "    $0 dir1 dir2 outdir"
    echo "  Combine output from dir1 and dir2 and put the aggregate"
    echo "  result in outdir.  Use this to achieve high phonon counts"
    echo "  by combining multiple smaller-count runs."
}

if [ $# -ne 3 ]; then
    echo "error: wrong number of arguments."
    print_usage;
    exit
elif [ ! -d "$1" ] || [ ! -d "$2" ]; then
    echo "error: dir1 and dir2 must be directories."
    print_usage;
    exit
elif [ -e "$3" ]; then
    echo "error: outdir must not exist."
    print_usage;
    exit
fi
if [ ! -f "$COMB_M" ]; then
    echo "error: can't find combine.m"
    exit
fi

# Create output dir:
dir1="$1"
dir2="$2"
odir="$3"
mkdir -p "$odir/seisfiles"

# Combo Log:
if [ -f "$dir1/combolog" ]; then
    cat "$dir1/combolog" >> "$odir/combolog"
else
    echo "$dir1" >> "$odir/combolog"
fi
if [ -f "$dir2/combolog" ]; then
    cat "$dir2/combolog" >> "$odir/combolog"
else
    echo "$dir2" >> "$odir/combolog"
fi

# Logfile combined:
if [ -f "$dir1/logfile-combined" ]; then
    cat "$dir1/logfile-combined" >> "$odir/logfile-combined"
else
    echo "$dir1" >> "$odir/logfile-combined"
    cat "$dir1/logfile" >> "$odir/logfile-combined"
    echo >> "$odir/logfile-combined"
fi
if [ -f "$dir2/logfile-combined" ]; then
    cat "$dir2/logfile-combined" >> "$odir/logfile-combined"
else
    echo "$dir2" >> "$odir/logfile-combined"
    cat "$dir2/logfile" >> "$odir/logfile-combined"
    echo >> "$odir/logfile-combined"
fi

# Combine:
for file1 in "$dir1"/$seispattern "$dir1"/$parampattern; do
    #bn=`basename $file1` # doesn't work if file is in subdirectory
    bn="${file1#$dir1}"
    file2="$dir2"/"$bn"
    ofile="$odir"/"$bn"
    if [ ! -f "$file2" ]; then
        echo "warning: files missing in dir2; skipping $bn."
        continue
    fi
    #echo "Combining: $file1 and $file2..."
    echo "combine(\"$file1\", \"$file2\", \"$ofile\")" \
        | octave -qfW -p "`dirname $COMB_M`"

done
