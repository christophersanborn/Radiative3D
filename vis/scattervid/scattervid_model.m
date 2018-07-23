## GNU/Octave script to make video frames from Radiative3D Scatter
## Data
##
## By Christopher J. Sanborn
##
##  This one plots from a side-elevation perspective, e.g, plot-y
##  corresponds to data-z, and plot-x corresponds to data-rho.  The
##  left side of the plot corresponds to the source epicenter and the
##  right side to some distant location, e.g., the location of a
##  seismic station.
##
##  Plots scatters against backdrop of model plot.  Model grid is
##  described in 'gridfile', and is assumed to be Cylinder type model
##  grid. Currently does not interpret WCGTetra grids (see
##  scattervid_p2p for an ad hoc script that does interpret WCGTetra.)
##
##  Destination point is specified as an [x y z] or as the filename of a
##  seisfile whose location is used as the destination point.
##
function scattervid_model (
             AnnoStruct = 0,        # Annotations struct (NOT USED YET)
             ArgStruct = struct()   # Additional optional args
           )
           # ArgStruct elements and defaults:
           #  destination - Loc of target station or seisfile containing same
           #  gridfile   - Name of file with grid data
           #  gridpair   - Which depth plumbs to use (ix values)
           #  duration   - Simtime to plot (lesser of this or PhononTTL)
           #  numframes  - Number of frames in video
           #  azifilt    - Only plot scatters within +/-Azi of SRC-DEST line
           #  mparams    - Metadata like source tensor
           #  seispat    - seis-file pattern (e.g. "seis_%03d.octv")
           #  seisrange  - which ones to plot (eg. [160,169:10:319])
           #  base       - Base location establishes the zero of the x axis,
           #               and base->dest it's range. Base defaults to source
           #               location, but can differ.
           #  window     - Viewport [left right top bot] Note: L and R are
           #               fractional and T and B are raw depths.
           #  txtkey     - Station names for annotation, e.g.:
           #               {{"LOP","MAK","WUS"}}. (Note: double braces required
           #               if struct() is used to create this field.
           #  loctext    - Name of source location. (Defaults to txtkey{1})

  destination = ProvideDefault(ArgStruct, "destination", [100 0 0]);
  gridfile = ProvideDefault(ArgStruct, "gridfile", "griddump.txt");
  gridpair = ProvideDefault(ArgStruct, "gridpair", [1 2]);
  Duration = ProvideDefault(ArgStruct, "duration", inf);
  NumFrames = ProvideDefault(ArgStruct, "numframes", 300); 
  mparamsfile = ProvideDefault(ArgStruct, "mparams", "out_mparams.octv");
  seispat = ProvideDefault(ArgStruct, "seispat", "");
  seisrange = ProvideDefault(ArgStruct, "seisrange", []);
  window = ProvideDefault(ArgStruct, "window", [-0.05,1.5,10,-130]);
  txtkey = ProvideDefault(ArgStruct, "txtkey", {"SRC","STA","STB"});
  LocText = ProvideDefault(ArgStruct, "loctext", txtkey{1});
  pixwidth = ProvideDefault(ArgStruct, "pixwidth", 1000);
  aspect = ProvideDefault(ArgStruct, "aspect", 15/9);
  AziFilt = ProvideDefault(ArgStruct, "azifilt", 25);
  
  # PARAMETERS:
  if (ischar(destination))
    SEIS = load(destination);
    if (isfield(SEIS,"SEIS")) SEIS=SEIS.SEIS; end # Catch old format
    DEST = SEIS.Location;     # Dest XYZ (Surface point)
  else
    DEST = destination;       # Else assume explicit destination
  end
  
  PP = load("scatters_P");
  SS = load("scatters_S");
  SRC = load("source_loc");
  BASE = [492.31000, -263.65000, 1.05000]; # "LOP" from the map
  META = load(mparamsfile);
  BASE = ProvideDefault(ArgStruct,"base", SRC);

  BASESRC = (SRC-BASE)(1:2);
  SRCDEST = (DEST-SRC)(1:2);
  SRCDESTHAT = SRCDEST / sqrt(sum(SRCDEST.*SRCDEST));
  BASEDEST = (DEST-BASE)(1:2);
  BASEDESTHAT = BASEDEST / sqrt(sum(BASEDEST.*BASEDEST));


  # Time Interval
  Duration = min(Duration, META.PhononTTL);
  dt = Duration/NumFrames;

  # Title
  ExEq = getexplfrac(META.EventSourceMT)
  EventText = "Event";
  if ExEq > 0.9
     EventText = "Explosion";
  elseif ExEq < 0.1
    EventText = "Earthquake";
  end
  Title = sprintf("%s at %s; Freq = %0.2f; %g sec", 
                  EventText, LocText, META.Frequency, Duration);


  # Process P arrays:

  P_T   = PP(:,1);
  P_X   = PP(:,2) - SRC(1);
  P_Y   = PP(:,3) - SRC(2);
  P_Z   = PP(:,4); # - SRC(3);
  P_XYZ = [P_X P_Y P_Z];
  P_XY  = [P_X P_Y];

  P_R   = sqrt(sum(P_XYZ.^2, 2));
  P_RHO = sqrt(sum(P_XY.^2, 2));
  P_PH  = atan2(P_Y, P_X);
  P_COS = (P_XY * SRCDESTHAT') ./ P_RHO; 

  p_x   = P_RHO .* sign(P_COS);
  p_y   = P_Z;

  P_TI  = floor(P_T./dt);               # Time Index

  # Process S arrays:

  S_T   = SS(:,1);
  S_X   = SS(:,2) - SRC(1);
  S_Y   = SS(:,3) - SRC(2);
  S_Z   = SS(:,4); # - SRC(3);
  S_XYZ = [S_X S_Y S_Z];
  S_XY  = [S_X S_Y];

  S_R   = sqrt(sum(S_XYZ.^2, 2));
  S_RHO = sqrt(sum(S_XY.^2, 2));
  S_PH  = atan2(S_Y, S_X);
  S_COS = (S_XY * SRCDESTHAT') ./ S_RHO; 

  s_x   = S_RHO .* sign(S_COS);
  s_y   = S_Z;

  S_TI  = floor(S_T./dt);               # Time Index


  ##
  ## Setup Paper Space:
  ##

  figwidth  = 20/3; # 6.666..., Reasonable for text size
  figheight = figwidth/aspect;
  paperdpi=pixwidth/figwidth;
  fonttitle  = 9.5;
  fontxylabel = 9.0;
  fontaxes = 9.0;

  figinit(figwidth,figheight,       # Clear and initialize fig in 
          "paperunits", "inches",   # remote headless-safe way
          "visible", "off");        #

  ##
  ## Setup Plot Space:
  ##

  plot_top = window(3);
  plot_bot = window(4);
  plot_interest_width = sqrt(sum((DEST(1:2)-BASE(1:2)).^2));
  plot_right = window(2) * plot_interest_width;
  plot_left = window(1) * plot_interest_width;


  ## Make baseplot with model on figure(1)

  set(gcf(), "visible", "off");
  hold off;
  plot (0,0); # dot at center, just to make axes command not get
              # clobbered
  axis([plot_left plot_right plot_bot plot_top]);
  set(gca(),"fontsize", fontaxes);
  box("off");
  title(Title, "fontsize", fonttitle);
  xlabel("Range (km)", "fontsize", fontxylabel);
  ylabel("Depth (km)", "fontsize", fontxylabel);
  hold on;
  modelbaseplot(gridfile,gridpair,mparamsfile,seispat,seisrange,txtkey);
  hBasePlot = gcf();  

  ##
  ## Make Frames:
  ##

  for f_idx = (0:(NumFrames-1))

    lstP = find( P_TI==f_idx & abs(P_COS)>=cos(AziFilt*(pi/180)) );
    lstS = find( S_TI==f_idx & abs(S_COS)>=cos(AziFilt*(pi/180)) );

    # Off-axis Correction

    OA_SX = (BASEDESTHAT*BASESRC') + s_x(lstS)*(BASEDESTHAT*SRCDESTHAT');
    OA_PX = (BASEDESTHAT*BASESRC') + p_x(lstP)*(BASEDESTHAT*SRCDESTHAT');
    OA_SY = s_y(lstS);  # This will eventually be un-flatten xform
    OA_PY = p_y(lstP);

    # Plot scatters:
    
    hFramePlot = copyobj(hBasePlot);    # Copy base plot to figure(2) and
    figure(hFramePlot);                 # select to overlay scatters
    set(gcf(), "visible", "off");       # Not sure why not copied from orig
    
    scatter(OA_SX, OA_SY, "b");
    scatter(OA_PX, OA_PY, "r");

    text(plot_right, plot_bot, sprintf("t = %5.2f s ",(f_idx+1)*dt),
              "verticalalignment", "bottom", "horizontalalignment", "right",
              "fontsize", fontxylabel, "fontweight", "demi");  # Time label

    filename = sprintf("framecache/FIG__%04d.png", f_idx);
    print(filename, sprintf("-r%1.0f",paperdpi));

    delete(hFramePlot);     # Clear room for next iteration

  endfor


  ##
  ## Produce Movie:
  ##

  file_cache = "framecache/FIG__%04d.png";
  file_out   = "framecache/movie.mp4";

  ffmakevid(file_out, file_cache, "pixwidth", pixwidth,
            "framerate", 10, "bitrate", 2000000); 

endfunction


function expl = getexplfrac(MT)

  gam = trace(MT);
  #lam = gam/3;
  #LAM = lam*eye(3);  ## The explosive (isotropic) part
  #SHR = MT - LAM;    ## The shear part

  explpart = (gam^2) / 3;   ## Mag^2 of explosive part
  total = sum(sum(MT.*MT)); ## Mag^2 of whole tensor

  expl = explpart/total;

end

function val = ProvideDefault(S, key, defvalue)
  if (isfield(S, key))
    val = getfield(S, key);
  else
    val = defvalue;
  end
end
