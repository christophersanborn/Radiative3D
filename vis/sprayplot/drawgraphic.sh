#!/bin/bash

DATANAME=points
PNGOUT=/scratch/dropbox/out_image.png

#PROJECTION='-Rg -JA20/25/6.5i'
PROJECTION='-R-90/270/-90/90 -JR6.5i'
ANOTATION='-B60g30/15g15 --GRID_PEN_PRIMARY=gray --ANNOT_FONT_SIZE_PRIMARY=9p'
SYMBOL=' -S -W4'
PAPEROPTS='-P'
DPI='160'

psxy $DATANAME $PROJECTION $ANOTATION $SYMBOL $PAPEROPTS -K > $DATANAME.ps 
echo "0 0" | \
psxy -O -K -R -J -Sc.15 -Ggreen >> $DATANAME.ps
echo "90 0" | \
psxy -O -K -R -J -Sc.15 -Gyellow >> $DATANAME.ps
echo "0 90" | \
psxy -O -R -J -Sc.15 -Gred >> $DATANAME.ps

ps2raster $DATANAME.ps -A -Tg -Qg4 -E$DPI -F$PNGOUT
