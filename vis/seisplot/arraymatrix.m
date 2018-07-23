## arraymatrix.m
##
## Takes an ARRAY struct containing a linear concatenation of seismic
## traces and returns a reshaped matrix suitable for plotting with
## image() or imagesc() plotting functions.
##
## optional arguments allow specification of which axes should be
## included in the output.
##
## Indexing the resulting matrix MM by MM(i,j), the i index indexes
## seismometers and the j index indexes time window.
##
function MM = arraymatrix(ARRAY, AXES=[1 1 1])

  MegaTrace = ARRAY.TraceXYZ * AXES';   # Results in a column vec

  SamplesPerSeis = length(ARRAY.TraceXYZ) / ARRAY.NumSeismometers;
  MM = reshape(MegaTrace, SamplesPerSeis, ARRAY.NumSeismometers)';

end
