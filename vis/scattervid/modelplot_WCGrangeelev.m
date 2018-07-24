# modelplot_WCGrangeelev.m
#
# Plots a range-elevation plot of the model, with annotations. Assumes
# grid is indexed in Range,Azimuth,Elevation, (i.e., that incrementing
# the first grid index takes you out a radial line. This is the
# structure used by the North Sea Crust Pinch models, for example.)
#
# Uses baseplot_gridWCGrangeelev.m for the actual plotting.  (The
# baseplot_ funciton is used both by this file and by
# scattervid_p2p.m)
#
function modelplot_WCGrangeelev(
             gridfile,            # Name of file with grid data
             param_reserved1 = 0, #
             plottitle = "Crust Model",     # Title
             mparfile = "out_mparams.octv", # MParams file
             seispattern = "",    # Glob pattern for seis files
             seisrange = [],      # Which ones to plot (eg. [160,169:10:319])
             PaperWH = [5.0 3.75],# Paper width, height
             wiremesh = true,     #
             RotMat = [1 0; 0 1]  # Optionals rotation/transform matrix, e.g.
                                  # to adjust orientation of an OCS RAW
                                  # "Earth-curved" grid
           )

  GG = read_gridgeom(gridfile);
  GG = GG(:,1,:,:);   # Discard azi index, not needed for range-elevation
  LayerColors = [ 0.5 0.3 0.1;
                  0.6 0.6 0.7;
                  0.7 0.1 0.1;
                  0.9 0.5 0.1; ];
  SEISMETA = read_seismeta(seispattern,seisrange);

  figinit(PaperWH(1), PaperWH(2),   # Initialize figure
          "paperunits", "inches",   #
          "visible", "off");        #

  rangemax = 0;
  depthmax = 0;
  [nR,nAz,nZ,nattr] = size(GG);
  for iR=1:nR
    for iZ=1:nZ
      R = sqrt(GG(iR,1,iZ,1)^2 + GG(iR,1,iZ,2)^2);
      Z = GG(iR,1,iZ,3);
      rangemax = max(R,rangemax);
      depthmax = min(Z,depthmax);
    end
  end
  viewwindow = [0,rangemax,depthmax-35,25];  # TODO: Smarten left and top

  hold on;
  plot (0,0);                     # Dot at center, (prevents
  axis(viewwindow);               #  clobbering of axis command).
  title(plottitle);               #
  xlabel("Range (km)");           #

  baseplot_gridWCGrangeelev(GG,LayerColors,RotMat); # Model color fill
  if (wiremesh)
    baseplot_gridWCGrangeelev(GG,0,RotMat);         # Plot (thicker) wire mesh
  end
  if (length(SEISMETA)>0)
    # TODO: Plot seismometers.  Note that scattervid_p2p defines an
    # inline function baseplot_seismap_elev to plot the seismometers.
    # This should probably be removed to a separate file and called
    # here.
  end


end
