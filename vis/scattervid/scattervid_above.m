## GNU/Octave script to make video frames from Radiative3D Scatter
## Data
##
## By Christopher J. Sanborn
##
##  This one plots from a plan-view or view-from-above perspective.
##  The viewport is computed to include the source as well as two
##  destination locations.
##

function scattervid_above(
             AnnoStruct = 0,  # Annotation struct (-1 signals lop nor defaults,
                              #  0 signals 'no annotations', otherwise see
                              #  ValidateAnnoStruct() for struct format.)
             ArgStruct = struct()       # Additional optional args
           )
           # ArgStruct elements and defaults:
           #  duration,  - Simtime to plot (lesser of this or PhononTTL)
           #  numframes, - Number of frames in video
           #  mparams,   - Metadata like source tensor
           #  gridfile   - Name of file with grid data
           #  loctext    - Name of src location; overrides AnnoStruct if set.
           #  seispat,   - seis-file pattern (e.g. "seis_%03d.octv")
           #  seisrange, - which ones to plot (eg. [160,169:10:319])
           #  window,    - Explicit viewport [X_west X_east Y_north Y_south]
           #               If set, overrides "smart" guess from GetViewExtents.

  Duration = ProvideDefault(ArgStruct, "duration", inf);
  NumFrames = ProvideDefault(ArgStruct, "numframes", 300); 
  mparamsfile = ProvideDefault(ArgStruct, "mparams", "out_mparams.octv");
  gridfile = ProvideDefault(ArgStruct, "gridfile", "griddump.txt");
  seispat = ProvideDefault(ArgStruct, "seispat", "");
  seisrange = ProvideDefault(ArgStruct, "seisrange", []);
  window = ProvideDefault(ArgStruct, "window", []);
  LocText = ProvideDefault(ArgStruct, "loctext", []);
  pixwidth = ProvideDefault(ArgStruct, "pixwidth", 1000);
  aspect = ProvideDefault(ArgStruct, "aspect", 4/3);
  

  PP = load("scatters_P");
  SS = load("scatters_S");
  SRC = load("source_loc");
  META = load(mparamsfile);
  SEISMETA = read_seismeta(seispat,seisrange);
  if (exist(gridfile,"file")==2)
    GG = read_gridgeom(gridfile);
    GG = GG(:,:,1,:);   # Discard depth index
    plotgrid = 1;       # 1: Plot WCG map view
    RangeColors = [ repmat([ 0.9 0.8 0.7 ], 5, 1);
                    repmat([ 0.7 0.5 0.4 ], 1, 1);
                    repmat([ 0.5 0.3 0.1 ], 2, 1);
                    repmat([ 0.7 0.5 0.4 ], 1, 1);
                    repmat([ 0.9 0.8 0.7 ], 3, 1);];
  else
    plotgrid = 0;       # 0: Plot nothing
  end
  AnnoStruct = ValidateAnnoStruct(AnnoStruct, SRC);

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
  if (ischar(LocText))
    EventLocText = LocText;
  elseif (length(AnnoStruct.SRCLabel)>0)
    EventLocText = AnnoStruct.SRCLabel;
  else
    EventLocText = "(Loc)";
  end
  Title = sprintf("%s at %s; Freq = %0.2f; %g sec", 
                  EventText, EventLocText, META.Frequency, Duration);


  # Process P arrays:

  P_T   = PP(:,1);
  P_X   = PP(:,2) ;# - SRC(1);
  P_Y   = PP(:,3) ;# - SRC(2);
  P_Z   = PP(:,4); # - SRC(3);
  P_XYZ = [P_X P_Y P_Z];
  P_XY  = [P_X P_Y];

  P_R   = sqrt(sum(P_XYZ.^2, 2));
  P_RHO = sqrt(sum(P_XY.^2, 2));
  P_PH  = atan2(P_Y, P_X);

  p_x   = P_X;
  p_y   = P_Y;

  P_TI  = floor(P_T./dt);               # Time Index

  # Process S arrays:

  S_T   = SS(:,1);
  S_X   = SS(:,2) ;# - SRC(1);
  S_Y   = SS(:,3) ;# - SRC(2);
  S_Z   = SS(:,4); # - SRC(3);
  S_XYZ = [S_X S_Y S_Z];
  S_XY  = [S_X S_Y];

  S_R   = sqrt(sum(S_XYZ.^2, 2));
  S_RHO = sqrt(sum(S_XY.^2, 2));
  S_PH  = atan2(S_Y, S_X);

  s_x   = S_X;
  s_y   = S_Y;

  S_TI  = floor(S_T./dt);               # Time Index


  ##
  ## Setup Paper Space:
  ##

  figwidth  = 6; # 
  figheight = figwidth/aspect;
  paperdpi=pixwidth/figwidth;
  #fonttitle  = 9.5;             # Not actually parameterized below, 
  fontxylabel = 9.0;             # but uncomment if I want to parameterize
  #fontaxes = 9.0;               # them later

  figinit(figwidth,figheight,       # Clear and initialize fig in 
          "paperunits", "inches",   # remote headless-safe way
          "visible", "off");        #


  # Compute some plot attributes:
  if (length(window)==4)
    viewwindow = window([1 2 4 3]); # [W E N S] to [L R B T]
  else
    viewwindow = GetViewExtents(AnnoStruct, GG);
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

  for f_idx = (0:(NumFrames-1))

    lstP = find(P_TI==f_idx);
    lstS = find(S_TI==f_idx);

    # Plot

    clf();
    hold on;
    if (plotgrid==1)                          #
      baseplot_gridWCGabove(GG,RangeColors);  # Plot grid mesh as backdrop
    end  
    if (length(SEISMETA)>0)             #
      baseplot_seismap_above(SEISMETA); # Plot seismometers (defined below)
    end
    if (plotgrid==1)                    #
      baseplot_gridWCGabove(GG);        # Plot grid mesh as backdrop
    end
    PlotAnnotations(AnnoStruct);        # Plot station labels and vectors
    scatter(s_x(lstS), s_y(lstS), "b");
    scatter(p_x(lstP), p_y(lstP), "r");
    title(Title);
    ylabel("Northing (km)");
    xlabel("Easting (km)");
    axis(viewwindow);

    text(plot_left+0.03*plot_width, plot_bot+0.97*plot_height,
         sprintf("t = %5.2f s ",(f_idx+1)*dt),
         "verticalalignment", "top", "horizontalalignment", "left",
         "fontsize", fontxylabel, "fontweight", "demi");  # Time label

    filename = sprintf("framecache/FIG__%04d.png", f_idx);
    print(filename, sprintf("-r%1.0f",paperdpi));

  endfor


  ##
  ## Produce Movie:
  ##

  file_cache = "framecache/FIG__%04d.png";
  file_out   = "framecache/movie.mp4";

  ffmakevid(file_out, file_cache, "pixwidth", pixwidth,
            "framerate", 10, "bitrate", 2000000); 


endfunction


########################################
### Helper Functions:
###

function RAS = ValidateAnnoStruct(AS, SRC)
  
  # First check for request-for-default:
  if (!isstruct(AS))  # If not a struct,
    Code = AS;        # then make note of what it is,
    AS = struct();    # and convert to struct.
    if (isscalar(Code) && (Code==-1))   # Then populate defaults for Lop Nor:
      AS.Dests      = [ -390.04 -167.18 1.457;      # Destination: WUS
                        -102.27  430.84 0.600  ];   # Destination: MAK
      AS.DestLabels = {"WUS", "MAK"};
      AS.AlignH     = {"center", "center"};
      AS.AlignV     = {"top", "bottom"};
      AS.SRCLoc     = SRC;
      AS.SRCLabel   = "Lop Nor";
      AS.Color      = [0 0.8 0.9];
    else  # Else populate with no-annotation defaults:
      AS.Dests = [];
      AS.SRCLoc = SRC;
      AS.SRCLabel   = "";
      AS.Color      = [0 0.8 0.9];
    end
  else
    if (!isfield(AS,"Dests")); AS.Dests = []; end
    if (!isfield(AS,"DestLabels")); AS.DestLabels = {}; end
    if (!isfield(AS,"AlignH")); AS.AlignH = {"center"}; end
    if (!isfield(AS,"AlignV")); AS.AlignH = {"top"}; end
    if (!isfield(AS,"SCRLoc")); AS.SRCLoc = SRC; end
    if (!isfield(AS,"SRCLabel")); AS.SRCLabel = {""}; end
    if (!isfield(AS,"Color")); AS.Color = [0 0.8 0.9]; end
  end

  RAS = AS;

end

function PlotAnnotations(AS)

  for iD = 1:(size(AS.Dests)(1))    # For each row in the Dests array
    X = [AS.SRCLoc(1), AS.Dests(iD,1)];
    Y = [AS.SRCLoc(2), AS.Dests(iD,2)];
    if (length(AS.DestLabels)>=iD)
      DLabel = AS.DestLabels{iD}; else DLabel = "";
    end
    if (length(AS.AlignH)>=iD)
      AlignH = AS.AlignH{iD}; else AlignH = "center";
    end
    if (length(AS.AlignV)>=iD)
      AlignV = AS.AlignV{iD}; else AlignV = "bottom";
    end
    line( X, Y, "LineWidth", 4, "Color", AS.Color);
    text(X(end), Y(end), DLabel,
         "color", AS.Color,
         "horizontalalignment", AlignH,
         "verticalalignment", AlignV);
  end
  X = [AS.SRCLoc(1)-1, AS.SRCLoc(1)+1]; # Small X at Source location
  Y = [AS.SRCLoc(2)-1, AS.SRCLoc(2)+1]; #
  line( X, Y, "LineWidth", 4, "Color", AS.Color);
  line( X, Y(end:-1:1), "LineWidth", 4, "Color", AS.Color);
  if (length(AS.SRCLabel) > 0)
       text(AS.SRCLoc(1), AS.SRCLoc(2), AS.SRCLabel,
            "color", AS.Color);
  end

end

function baseplot_seismap_above(SM)

  for idx = 1:length(SM)
      GatherAx1 = SM(idx).GatherRadius(1,2) * SM(idx).AxesX1;
      GatherAx2 = SM(idx).GatherRadius(1,2) * SM(idx).AxesX2;
      LineColor = [0.1 0.4 0.6];
      FillColor = [0.6 0.9 0.6];
      LineWidth = 1;
      TwoAxEllipse(GatherAx1, GatherAx2, SM(idx).Location,
                   LineWidth, FillColor, FillColor);
      X = SM(idx).Location(1) + [-0.7*GatherAx1(1) -GatherAx1(1)];
      Y = SM(idx).Location(2) + [-0.7*GatherAx1(2) -GatherAx1(2)];
      line (X, Y, "linewidth", LineWidth,
                  "color", LineColor);
      X = SM(idx).Location(1) + [-0.7*GatherAx2(1) -GatherAx2(1)];
      Y = SM(idx).Location(2) + [-0.7*GatherAx2(2) -GatherAx2(2)];
      line (X, Y, "linewidth", LineWidth,
                  "color", LineColor);
      X = SM(idx).Location(1) + [-0.2*GatherAx1(1) GatherAx1(1)];
      Y = SM(idx).Location(2) + [-0.2*GatherAx1(2) GatherAx1(2)];
      line (X, Y, "linewidth", LineWidth,
                  "color", LineColor);
      X = SM(idx).Location(1) + [-0.2*GatherAx2(1) GatherAx2(1)];
      Y = SM(idx).Location(2) + [-0.2*GatherAx2(2) GatherAx2(2)];
      line (X, Y, "linewidth", LineWidth,
                  "color", LineColor);
  end

end

function TwoAxEllipse(V1, V2, Center, LineWidth, LineColor, FillColor)

  persistent ucosphi; # Unit cosines over 2pi range of phi's
  persistent usinphi; #  ''   sines
  persistent binit = false;

  if (!binit)         # Pre-compute a bunch of trig calculations.
    phis=linspace(0,2*pi,33); # (not sure if major time saver)
    ucosphi=cos(phis);
    usinphi=sin(phis);
    binit=true;
  end

  X = V1(1)*ucosphi + V2(1)*usinphi + Center(1);
  Y = V1(2)*ucosphi + V2(2)*usinphi + Center(2);

  h = fill (X, Y, FillColor);
  set (h, "linewidth", LineWidth, 
          "edgecolor", LineColor
       );

end

function VE = GetViewExtents(AS, GG, Aspect=1.3, MarginFactor=1.05)

  if (ndims(GG)==4) # Then get view extents from grid
    plot_inner_bound_maxX = max(max(GG(:,:,1,1)));
    plot_inner_bound_minX = min(min(GG(:,:,1,1)));
    plot_inner_bound_maxY = max(max(GG(:,:,1,2)));
    plot_inner_bound_minY = min(min(GG(:,:,1,2)));     
  else              # Else get it from DESTS
    plot_inner_bound_maxX = max([AS.SRCLoc(1); AS.Dests(:,1)]);
    plot_inner_bound_minX = min([AS.SRCLoc(1); AS.Dests(:,1)]);
    plot_inner_bound_maxY = max([AS.SRCLoc(2); AS.Dests(:,2)]);
    plot_inner_bound_minY = min([AS.SRCLoc(2); AS.Dests(:,2)]);
  end 
  plot_midX = (plot_inner_bound_minX + plot_inner_bound_maxX)/2; 
  plot_midY = (plot_inner_bound_minY + plot_inner_bound_maxY)/2; 
  plot_inner_widX = plot_inner_bound_maxX - plot_inner_bound_minX;
  plot_inner_widY = plot_inner_bound_maxY - plot_inner_bound_minY;
  if (plot_inner_widX < Aspect*plot_inner_widY)
     plot_inner_widX = Aspect*plot_inner_widY;
  else
     plot_inner_widY = plot_inner_widX/Aspect;
  end
  plot_widX = MarginFactor * plot_inner_widX;
  plot_widY = MarginFactor * plot_inner_widY;
  plot_left = plot_midX - 0.5*plot_widX;
  plot_right = plot_midX + 0.5*plot_widX;
  plot_top = plot_midY + 0.5*plot_widY;
  plot_bot = plot_midY - 0.5*plot_widY;

  VE = [plot_left plot_right plot_bot plot_top];

end

function expl = getexplfrac(MT)

  gam = trace(MT);
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
