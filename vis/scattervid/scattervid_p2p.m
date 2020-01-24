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
##  NOTE: This should probably be renamed to scattervid_wcgrae or
##  something, as this basically does the same a scattervid_model but
##  for WCGTetra grids arranged as range, azi, elev, instead of for
##  cylinder grids.  There are some other difference, like not
##  providing for the elevation plane to differ from the SRC->DEST
##  path. (c.f. _model's handling of BASE->DEST vs SRC->DEST.)  But
##  these scripts are otherwise now very similar.
##
##  NOTE:  At the moment, this script does NOT filter by azimuth.  It
##  takes the whole azimuth and compresses into rho.  But it is
##  *intended* to filter by azimuth, eventually.
##
##  If called with non-empty SelectedIdx, can be used to print
##  individual frames as pdf's, instead of producing a whole movie.
##
##  Return values are handles to the scatter point series (of the most
##  recent frame).  Can be used to tweak how they display (when used
##  to produce a single frame).
##
function [hP hS] = scattervid_p2p(
             AnnoStruct = 0,        # Annotations struct (NOT USED YET)
             ArgStruct = struct(),  # Additional optional args
             PaperWH = [20/3],      # [Width, Height]. Height is assumed from
                                    # aspect ratio if absent.
             SelectedIdx = []       # If empty we produce movie. If a list, then
                                    # no movie but we print individual frames as
                                    # pdfs. If length one we plot but do not
                                    # print (assumed top level will tweak titles
                                    # and such before print).
           )
           # ArgStruct elements and defaults:
           #  destination - Loc of target station or seisfile containing same
           #  gridfile   - Name of file with grid data
           #  gridcolors - Layer colors as column of RGB row-triples
           #  pmarkprops -
           #  smarkprops -
           #  axesprops  -
           #  duration   - Simtime to plot (lesser of this or PhononTTL)
           #  numframes  - Number of frames in video
           #  pixwidth   -
           #  azifilt    - Only plot scatters within +/-Azi of SRC-DEST line
           #  mparams    - Metadata like source tensor
           #  seispat    - seis-file pattern (e.g. "seis_%03d.octv")
           #  seisrange  - which ones to plot (eg. [160,169:10:319])
           #  loctext    - Source location name
           #  titlefmt   - Title fprint string; Overrides loctext if set
           #  window     - Explicit viewport [L R T B]
           #               If set, overrides "smart" best guess
           #  fileprefix - Prefix on frame files when individual frames requested
           #
           # Note on marker props: Pass this as a double-braced array. E.g.
           # Args.pmarkprops = {{"linewidth", 2.5[, ...]}}.  Also, props are
           # for a scatter series, not a plot line. So use "sizedata" instead
           # of "size", and there is not "color", etc.

  destination = ProvideDefault(ArgStruct, "destination", [100 0 0]);
  gridfile = ProvideDefault(ArgStruct, "gridfile", "griddump.txt");
  Duration = ProvideDefault(ArgStruct, "duration", inf);
  NumFrames = ProvideDefault(ArgStruct, "numframes", 300); 
  mparamsfile = ProvideDefault(ArgStruct, "mparams", "out_mparams.octv");
  seispat = ProvideDefault(ArgStruct, "seispat", "");
  seisrange = ProvideDefault(ArgStruct, "seisrange", []);
  window = ProvideDefault(ArgStruct, "window", []);
  LocText = ProvideDefault(ArgStruct, "loctext", "(Loc)");
  TitleFTxt = ProvideDefault(ArgStruct, "titlefmt", []);
  pixwidth = ProvideDefault(ArgStruct, "pixwidth", 1000);
  aspect = ProvideDefault(ArgStruct, "aspect", 5/3);
  LayerColors = ProvideDefault(ArgStruct, "gridcolors", []);
  pmarkprops = ProvideDefault(ArgStruct, "pmarkprops", {});
  smarkprops = ProvideDefault(ArgStruct, "smarkprops", {});
  axesprops = ProvideDefault(ArgStruct, "axesprops", {});
  AziFilt = ProvideDefault(ArgStruct, "azifilt", 25); ## Not Used
  fileprefix = ProvideDefault(ArgStruct, "fileprefix", "");

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
  META = load(mparamsfile);
  SEISMETA = read_seismeta(seispat,seisrange);
  if (exist(gridfile,"file")==2)
    GG = read_gridgeom(gridfile);
    GG = GG(:,1,:,:);   # Discard azi index, not needed for range-elevation
    plotgrid = 1;       # 1: Plot WCG range-elevation profile
  else
    plotgrid = 0;       # 0: Plot nothing
  end
  if (size(LayerColors,2)!=3)     # If not an Nx3 matrix...
    LayerColors = [ 0.5 0.3 0.1;  # Then provide a default value.
                    0.6 0.6 0.7;
                    0.7 0.1 0.1;
                    0.9 0.5 0.1;];
  end

  # Time Interval
  Duration = min(Duration,META.PhononTTL);
  dt = Duration/NumFrames;

  # Title
  ExEq = getexplfrac(META.EventSourceMT)
  EventText = "Event";
  if ExEq > 0.9
     EventText = "Explosion";
  elseif ExEq < 0.1
    EventText = "Earthquake";
  end
  Title = sprintf("%s at %s; Freq = %0.3f; %g sec", 
                  EventText, LocText, META.Frequency, Duration);
  if (ischar(TitleFTxt))
    Title = sprintf(TitleFTxt, META.Frequency, Duration);
  end

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

  p_x   = P_RHO;
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

  s_x   = S_RHO;
  s_y   = S_Z;

  S_TI  = floor(S_T./dt);               # Time Index


  ##
  ## Setup Paper Space:
  ##

  figwidth  = PaperWH(1);           # (Default 6.666_, reasonable for text size)
  figheight = figwidth/aspect;      # Height automatic,
  if (length(PaperWH>1))            # unless explicit.
    figheight = PaperWH(2);
  end
  paperdpi=pixwidth/figwidth;
  #fonttitle  = 9.5;    # TODO: font sizes not parameterized yet
  fontxylabel = 9.0;
  #fontaxes = 9.0;

  figinit(figwidth,figheight,       # Clear and initialize fig in
          "paperunits", "inches",   # remote headless-safe way
          "visible", "off");        #

  ##
  ## Setup Plot Space:
  ##

  plot_top = max(max(P_Z),max(S_Z)) + 25;
  plot_bot = min(min(P_Z),min(S_Z)) - 35;
  plot_width = 1.7 * sqrt(sum((DEST(1:2)-SRC(1:2)).^2));
  viewwindow = [0 plot_width plot_bot plot_top];
  if (length(window)==4)   # If user passed a 'window' arg then override
    viewwindow = window([1 2 4 3]);     # [LRTB] to [LRBT]
  end
  plot_left = viewwindow(1);
  plot_right = viewwindow(2);
  plot_top = viewwindow(4);
  plot_bot = viewwindow(3);
  plot_width = plot_right-plot_left;
  plot_height = plot_top-plot_bot;

  ##
  ## Make Frames:
  ##

  makemovie_b = true;                   # << Behavior flags
  printframes_b = true;                 #
  f_idx_range = (0:(NumFrames-1));
  if (length(SelectedIdx)>0)
    f_idx_range = SelectedIdx;
    makemovie_b = false;
    if (length(SelectedIdx)==1) printframes_b = false; end
  end

  for f_idx = f_idx_range               # LOOP over frames:

    lstP = find(P_TI==f_idx);
    lstS = find(S_TI==f_idx);

    # Plot

    set(gcf(), "visible", "off");       # Initialize figure
    clf(); hold on;                     # Clear and hold
    plot (0,0);                         # Dot at center, (prevents
    axis(viewwindow);                   # clobbering of axes command)
    title(Title);                       # Titles and labels
    xlabel("Range (km)");               #

    if (plotgrid==1)
      baseplot_gridWCGrangeelev(GG,LayerColors);  # Model color fill
      baseplot_gridWCGrangeelev(GG);              # Plot grid wire mesh
    end
    if (length(SEISMETA)>0)
      baseplot_seismap_elev(SEISMETA);            # Plot seismometers (defined below)
    end

    hS2 = scatter(s_x(lstS), s_y(lstS), 8, "none"); # Plot Phonons: shadow
    hS = scatter(s_x(lstS), s_y(lstS), 4, "none");  #               marker
    hP2 = scatter(p_x(lstP), p_y(lstP), 8, "none"); # (TODO: would 'plot' be
    hP = scatter(p_x(lstP), p_y(lstP), 4, "none");  # (better than 'scatter'?)
    set(hP2,"linewidth",1.0,"markerfacecolor",[0.4 0 0],"marker","s");
    set(hP,"linewidth",1.0,"markerfacecolor",[1.0 0 0],"marker","s");
    set(hS2,"linewidth",1.0,"markerfacecolor",[0 0 0.4],"marker","s");
    set(hS,"linewidth",1.0,"markerfacecolor",[0 0 1.0],"marker","s");
                                                    # Sqr mrkr reduces filesize
    text(plot_right, plot_bot,
              sprintf("t = %5.2f s ",(f_idx+1)*dt), "tag", "TimeCode",
              "verticalalignment", "bottom", "horizontalalignment", "right",
              "fontsize", fontxylabel, "fontweight", "normal");  # Time label

    set(hP2, pmarkprops(1:2:end), pmarkprops(2:2:end));
    set(hP, pmarkprops(1:2:end), pmarkprops(2:2:end));
    set(hS2, smarkprops(1:2:end), smarkprops(2:2:end));
    set(hS, smarkprops(1:2:end), smarkprops(2:2:end));
    set(gca(), axesprops(1:2:end), axesprops(2:2:end));
    if (makemovie_b)
      filename = sprintf("framecache/FIG__%04d.png", f_idx);
      print(filename, sprintf("-r%1.0f",paperdpi));
    else
      if (printframes_b)
        filename = sprintf("%sScatvid_P2P__Frame_%04d.pdf", fileprefix, f_idx);
        print(filename);
      end
    end

  endfor

  ##
  ## Produce Movie:
  ##

  if (makemovie_b)
    file_cache = "framecache/FIG__%04d.png";
    file_out   = "framecache/movie.mp4";

    ffmakevid(file_out, file_cache, "pixwidth", pixwidth,
              "framerate", 10, "bitrate", 2000000); 
  end

endfunction


############################################
### Helper Functions:   (scattervid_p2p.m)
###

function baseplot_seismap_elev(SM)

  for idx = 1:length(SM)

      StickLength = SM(idx).GatherRadius(1,2);
      StickLength = 6; # Override
      RZLoc = RhoZ(SM(idx).Location);
      RZAx3  = RhoZ(SM(idx).Location 
                    + 1.0*StickLength*SM(idx).AxesX3);
      RZAx3N = RhoZ(SM(idx).Location 
                    - 0.2*StickLength*SM(idx).AxesX3);
      X = [RZAx3N(1) RZAx3(1)];
      Y = [RZAx3N(2) RZAx3(2)];
      line(X, Y, "linewidth", 2, "color", [0.1 0.6 0.2]);
      scatter(RZLoc(1), RZLoc(2), "green");

  end

end

function rz = RhoZ(V)
  rz = [sqrt(sum(V(1:2).^2)),V(3)];
end

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
