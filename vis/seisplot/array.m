## array.m
##
## Reads a range of Seismometer files (containing structures defining
## metadata and trace data for a seismometer output from a Radiative3D
## run) and joins them together into a structure that defines a
## travel-time curve array.  This structure contains, essentially, a
## concatenation of the trace data, and metadata about the array, such
## as the number of seismometers included, min and max distances, time
## window, azimuth, etc.
##
## Takes a filename pattern, in printf format, and two integers
## denoting the beginning and end of a range of filenames.  Outputs an
## ARRAY struct.
##
##
##   Example:  AR = array("seis_%03d.octv", 0, 7)  
##
##
function ARRAY = array(fpat, ibegin, iend)

  ARRAY.NumSeismometers = 0;    # (increments below)
  ARRAY.TimeWindow = [];
  ARRAY.Distances = [];
  ARRAY.Azimuths = [];
  ARRAY.TraceXYZ = [];

  for i = ibegin:iend

    ARRAY.NumSeismometers += 1;
    load(sprintf(fpat,i));      # Loads a struct called 'SEIS'

    ARRAY.TraceXYZ = [ARRAY.TraceXYZ; SEIS.TraceXYZ];
    dist = range_km(SEIS.Location, SEIS.EventLoc);
    ARRAY.Distances = [ARRAY.Distances; dist];
    azimuth = azimuth_deg(SEIS.Location, SEIS.EventLoc);
    ARRAY.Azimuths = [ARRAY.Azimuths; azimuth];
    if i==ibegin
      ARRAY.TimeWindow = SEIS.TimeWindow;
      ARRAY.NumBins = SEIS.NumBins;
      ARRAY.Frequency = SEIS.Frequency;
      ARRAY.EventMT = SEIS.EventMT;
    end

  end


end
