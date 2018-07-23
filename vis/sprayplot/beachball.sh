#!/bin/bash
##
## Plots the P-phonon density using roughly the same viewpoint and
## projection scheme as the standard "beach ball" plot. 
##
## The standard beach ball plot places a sphere around the source, and
## marks P-polarity (+/-) on the sphere.  The plot conceptually
## consists of removing the upper hemisphere and viewing the bottom
## hemisphere from above.
##
## Our plot marks P-energy magnitude in the same manner.  (So we lose
## the polarity information, but we can live without it.)
##
## The data, before preprocessing, records energy output in "longitude
## and lattitude" about a sphere surrounding the source.  The sphere
## uses the coordinate system used by Aki and Richards, which has the
## +Z axis pointing downwards, the +X pointing north, and the +Y
## pointing east.
##
## Because our plotting algorithm can't plot the "inside" of the
## sphere, we have to first reflect the data during preprocessing and
## then plot a view of the outside of the sphere.  Because Aki's
## system points down, we actually plot the upper hemisphere of our
## sphere.  We orient the sphere with +Z upwards (out of page), +X
## northwards (top of page), and +Y westwards (left of page).  In
## preprocesing, we negate the azimuthal angle, which effectively
## inverts east and west, which gives us the view that we want (since
## +Y is SUPPOSED to point eastwards).
##
## Thus, when you look at the plot, you actually see NESW pointing Up,
## Right, Down, Left, as expected, and the density represent the
## bottom hemisphere surrounding the source, as it should.
##

DATANAME=points_refl   ## Preprocessed - reflected and S removed 
PNGOUT=beachball.png

PROJECTION='-Rg -JS180/89.99/6.5i'
ANOTATION='-B60g30/15g15 --GRID_PEN_PRIMARY=gray --ANNOT_FONT_SIZE_PRIMARY=9p'
SYMBOL=' -S -W4'
PAPEROPTS='-P'
DPI='160'

psxy $DATANAME $PROJECTION $ANOTATION $SYMBOL $PAPEROPTS `# -K` > $DATANAME.ps 

## Plot dots to help clarify axis directions:
## (Uncomment if wanted for diagnostic purposes and also add -K option to previous command above)
#
#echo "0 1" | \
#psxy -O -K -R -J -Sc.15 -Ggreen >> $DATANAME.ps     # North (Green / +X-axis)
#echo "90 1" | \
#psxy -O -K -R -J -Sc.15 -Gyellow >> $DATANAME.ps    # East (Yellow / +Y-axis)
#echo "0 90" | \
#psxy -O -R -J -Sc.15 -Gred >> $DATANAME.ps          # Down (Red / +Z-axis)
#

ps2raster $DATANAME.ps -A -Tg -Qg4 -E$DPI -F$PNGOUT
