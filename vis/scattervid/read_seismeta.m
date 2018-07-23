# read_seismeta.m
#
# Reads Seismometer file metadata (i.e. everything but the traces)
# into a Struct Array.  This can be used, e.g., to plot the
# seismometers on a map.
#
# The output is a 1xN matrix of Structs containing the fields (xxx).
# Access is by index and filedname, e.g.:  SM(n).Location
#
function SM = read_seismeta(
                seispat = "",   # seis-file pattern (e.g. "seis_%03d.octv")
                seisrange = []  # which ones to plot (eg. [160,169:10:319])
                )

  idx = 0;
  SM = [];
  for j=seisrange
    seisname=sprintf(seispat,j);
    if (exist(seisname,"file")==2)
      idx++;
      SEIS=load(seisname);
      if (isfield(SEIS,"SEIS")) SEIS=SEIS.SEIS; end # Catch old format
      SM(idx).Location = SEIS.Location;
      SM(idx).AxesX1 = SEIS.AxesX1;
      SM(idx).AxesX2 = SEIS.AxesX2;
      SM(idx).AxesX3 = SEIS.AxesX3;
      SM(idx).AxesDesc = SEIS.AxesDesc;
      SM(idx).GatherRadius = SEIS.GatherRadius;
    end
  end


end
