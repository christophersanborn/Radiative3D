# seisplot.m
#
# Makes an envelope plot from the trace data in 'tracefile', plotting
# X,Y,Z, and P, and S components.  Optionally applies gamma scaling to
# act as a sort of contrast control.  The 'gpow' argument determines
# this scaling, with:
#
#   gpow = 1.0   -->   Envelopes proportional to amplitude
#   gpow = 2.0   -->   Envelopes proportional to energy
#   gpow = X     -->   Envelop proportional to trace data^(X/2).
#                      Low values bring out detail in the lows.  High
#                      values enhance the peaks.
#
# multiseis.m
#
#
#
# Printing to file:
#
#   print("fig.png", "-r160");  # for  800x600 png file
#   print("fig.png", "-r240");  # for 1200x900 png file
#
#   Do not use the "-S" option to print - it won't do what you want.
#
# Note: Plot visibility is set to "off" by default to speed up
# scripted runs.  If making plot interactiviely, may need to
# set(gcf(),"visible", "on") after this function in order to see it.
#
# PARAMETERS:
#
#   vports = {[height, above, below]}
#
#     Defines the vertical height requirements of each trace.
#     vports{1} is for the topmost trace, vports{n} for the bottommost
#     trace. 'height' is the height of teh positive portion of the Y
#     axis for the trace.  'above' and 'below' define the minimum
#     needed buffer space needed above or below the positive-Y axis.
#     The positive range of any traces will be scaled to fit into
#     'height' space, and any negative values will swing below, so it
#     is important to ensure that 'below' is big enough to
#     accommodate. The gap between {n} and {n+1} will be
#     max(vports{n}(3),vports{n+1}(2).
#
#   vpwidth = [width, tmin, tmax]
#
#     Whatever 'width' is, the viewport will still span the width of
#     the figure, but it will serve to establish the relative meaning
#     of 'height' in the vports{} param.  'tmin' and 'tmax' establish
#     the common X-axis range for all viewports, typically in seconds.
#
#   tracefiles = {filename}
#
#     A cell array containing filenames of seismic trace files. Files
#     in the seisfile format are supported.  (Each such file contains
#     7 traces: X Y Z, P S, and CP CS.)  A simple one-trace format
#     will also be supported.
#
#   traces = {[tindex vpindex gpow ymax color blcolor fillcolor thick blthick]}
#
#     'tindex': index into traces found in the tracefiles.  'vpindex':
#     which viewport to plot in.  color args: index into colortable, 0
#     means don't show, -1 means default.
#
function multiseis (tracefiles,
                    vports={[1 0.8 .05],
                            [1 .05 .05],
                            [1 .05 .05],
                            [1 .05 .05],
                            [1 .05 1.1]},
                    vpwidth = [12 0 600],
                    traces = {[6 4  2.0 0  0 0 6  3.75 1.5],
                              [7 5  2.0 0  0 0 7  3.75 1.5],
                              [1 1  1.0 0  1 1 0  3.75 1.5],
                              [2 2  1.0 0  2 2 0  3.75 1.5],
                              [3 3  1.0 0  3 3 0  3.75 1.5],
                              [4 4  1.0 0  4 4 0  3.75 1.5],
                              [5 5  1.0 0  5 5 0  3.75 1.5],},
                    annotations = {},
                    PlotLabels = struct("Title","[Title]")
                  )

  ## Globals:
  global VP;            # Viewport parameters
  VP.vports = vports;   #
  VP.vpwidth = vpwidth; #

  global TR;            # Trace list (cell array of structs)
  GetTraces(tracefiles) # Populate


  ##
  ## Plot Prep: Layout and Style
  ##

  fonttitle  = 9.5;
  fontxylabel = 9;
  fontaxes  = 9;
  
  axwidth = vpwidth(1);             # Axes relative width
  axheight = GetTotalVPHeight();    # Also initializes various VP arrays
  axaspect = axheight/axwidth;      # H/W aspect (NOT the standard W/H)
  figheight = 4.5*axaspect+1;       # Fig H and W, assume +1" vertical for 
  figwidth = 5;                     # titles and labels and +0.5" horiz.

  figinit(figwidth, figheight,      # Clear and initialize 
       "paperunits", "inches",      # a 5.0" x 3.75" figure
       "visible", "off");
  axes("fontweight", "bold",        # Create and setup an axes object
       "fontsize", fontaxes,        #
       "linewidth", 1.5, "ytick", VP.lstY0,
       "yticklabel", {},
       "nextplot", "add");  # equiv to "hold on"

  axis([VP.vpwidth(2) VP.vpwidth(3) 0 VP.Ytotal]);

  ## 
  ## Begin Plotting:
  ##

  ## Iterate over trace array
  for i = 1:length(traces)

       Trace = TR{traces{i}(1)};
       ivp = traces{i}(2);
       gpow = traces{i}(3);     # Desired scale (1.0: amp 2.0: enrgy)
       Tgpow = Trace.Tgpow;     # Scale trace was encoded with
       YFullScale = traces{i}(4);
       color     = ColorMap(traces{i}(5));
       blcolor   = ColorMap(traces{i}(6));
       fillcolor = ColorMap(traces{i}(7));
       lthick  = traces{i}(8);
       blthick = traces{i}(9);

       if (YFullScale==0)
         YFullScale = Trace.ChFamMax ^ (gpow/Tgpow);
       elseif (YFullScale<0)
         YFullScale = (-YFullScale)*(Trace.ChanMax^(gpow/Tgpow));
       end # Else YFullScale is explicit; gpow assumed applied.

       Data = Trace.Data;
       Data = sign(Data) .* (abs(Data) .^ (gpow/Tgpow));

       VPplot(Trace.XRange, Data, ivp,
              color, blcolor, fillcolor,
              lthick, blthick, YFullScale);
  end


  ##
  ## Text Annotations:
  ##

  # Titles and axes labels:
  if (isfield(PlotLabels,"Title"))
    title(PlotLabels.Title, "fontsize", fonttitle);
  end
  xlabel("Time (s)", "fontsize", fontxylabel);
  if (isfield(PlotLabels,"YTickLabel"))
     set(gca(),"yticklabel",PlotLabels.YTickLabel);
  end

  for i=1:length(annotations)
      Anno = annotations{i};
      if(!isfield(Anno,"ivp"))  # PITA error checking here:
        error("Annotation structure MUST contain field 'ivp'.");
      end                       # Test that ivp filed exists,
      try                       # AND
        errmsg="";              #  ...
        if(!isscalar(Anno.ivp)) # that it is scalar,
          errmsg="Annotation ivp field must be a scalar integer.";
        end                     # AND catch the "unexpected behavior" of the 
      catch                     # struct() function when passed a cell array.
        error(sprintf("%s%s%s", # (Pass {cell_arr} or {{...}} instead.)
          "Whoops! Annotation struct got vectorized.\n",
          "\tYou likely passed a cell-array to the struct() function.\n",
          "\tTry again with an extra layer of curly braces {}."));
      end_try_catch             # Effin bloody octave hell...
      if (length(errmsg)>0)     # Can't call error in the try block or it
         error(errmsg);         # "catches" my own error... so we keep msg 
      end                       # for after the try block. Srsly f'n bldyoctvel
      msg = ""; if(isfield(Anno,"Text")); msg = Anno.Text; end;
      ivp = 1;  if(isfield(Anno,"ivp")); ivp = Anno.ivp; end;
      dx = 0;   if(isfield(Anno,"dX")); dx = Anno.dX; end;
      dy = 0;   if(isfield(Anno,"dY"));  dy = Anno.dY; end;
      side = "left"; if(isfield(Anno,"Side")); side = Anno.Side; end;
      pos = "base";  if(isfield(Anno,"Pos")); pos = Anno.Pos; end;
      tsize = 7; if(isfield(Anno,"Size")); tsize = Anno.Size; end;
      if(isfield(Anno,"Props")); tsize = Anno.Props; end; # if Props set, it
                                                            # must specify size
      if(length(msg)>0)
        VPannotate(msg,ivp,dx,dy,side,pos,tsize);
      end
      if(isfield(Anno,"LineXYXY"))
        lspec = {}; 
        if (isfield(Anno,"LProps")); lspec = Anno.LProps; end
        VPline(Anno.LineXYXY, ivp, lspec);
      end
  end

end# END main function
###
##

##**************************** 
##**  ===================
##**   HELPER FUNCTIONS:
##**  ===================
##*

##
##
function h = fillplot(X,Y,C,Y0=0)

  X = reshape(X,length(X),1);   # Ensure column vectors
  Y = reshape(Y,length(Y),1);   #

  # Fixes render bug when clipped extremities are non zero:
  Xmax=axis()(2);  # Only correct if intended axis range already set...
  i=find(X<=Xmax,1,"last");
  X=X(1:i);        # Truncate series just at or before edge of window.
  Y=Y(1:i);        # 
  Xmin=axis()(1);  # And repeat, but this time truncate left edge:
  i=find(X>=Xmin,1,"first");
  X=X(i:end);
  Y=Y(i:end);

  h = fill([X(1); X; X(end)],
           [  Y0; Y; Y0    ], C);
  set(h, "linestyle","none", "edgecolor", C);

end


######
##
function expl = getexplfrac(MT)

  tr = trace(MT);
  #lam = gam/3;
  #LAM = lam*eye(3);  ## The explosive (isotropic) part
  #SHR = MT - LAM;    ## The shear part

  explpart = (tr^2) / 3;   ## Mag^2 of explosive part
  total = sum(sum(MT.*MT)); ## Mag^2 of whole tensor

  expl = explpart/total;
  if tr<0; expl*=-1; end

end


######
## Calculates total Y height of all viewports and creates arrays of Y0
## values, etc.
function height = GetTotalVPHeight()
  global VP;
  tmpvp = {[0,0,0], VP.vports{end:-1:1}};   # Prepend dummy to rev'd vport list
  Ytotal=0;                             # Accumulator for total Y
  Y0rlist=[];                           # Bottom to top Y0 list (reversed)
  Y1rlist=[];                           # Bot to top Y1 list, where Y1=Y0+dY
  dYrlist=[];                           #
  for i = 1:length(VP.vports)
    j = i+1;
    belowskip = max(tmpvp{j}(3),tmpvp{i}(2));
    portskip = tmpvp{j}(1);
    Ytotal+=belowskip;
    Y0rlist(end+1)=Ytotal;
    dYrlist(end+1)=portskip;
    Ytotal+=portskip;
    Y1rlist(end+1)=Ytotal;
  end
  Ytotal+=tmpvp{end}(2);
  VP.Ytotal = Ytotal;
  VP.lstY0 = Y0rlist(end:-1:1); # Use lstXX for top-to-bot lists
  VP.lstDY = dYrlist(end:-1:1); #
  VP.lstY1 = Y1rlist(end:-1:1);
  VP.rlstY0 = Y0rlist;          # Use rlstXX for bot-to-top lists
  VP.rlstDY = dYrlist;          #
  VP.rlstY1 = Y1rlist;
  height = Ytotal;
end


######
##
function VPplot(Xdat, Ydat, ivp, lrgb, blrgb, fillrgb, lthk, blthk, Ymax=0)

  global VP;

  Y0 = VP.lstY0(ivp);
  DY = VP.lstDY(ivp);

  if (length(Xdat)==2)  # Allow Xdat = [min max], else expect it to match
     Xdat = linspace(Xdat(1),Xdat(2),length(Ydat));   # dimension of Ydat
  end                   

  # Ymax: (The Ydat value that maps to Y0+DY)
  if (Ymax==0); Ymax = -1; end  # Zero -> "default"
  if (Ymax<0);                  # Negative implies relative to 
    Ymax = (-Ymax)*max(Ydat);   #   found maximum.
  end                           # Positive imlpies absolute

  # Scale Ydat:
  Ydat = DY*Ydat./Ymax;

  # Fill, if requested:
  if (!isscalar(fillrgb))
    fillplot(Xdat,Y0+Ydat,fillrgb,Y0);
  end

  # Baseline, if requested:
  if (!isscalar(blrgb))
    FlatlineX = [VP.vpwidth(2) VP.vpwidth(3)];
    plot(FlatlineX,Y0+[0,0],"linewidth",blthk,"color",blrgb);
  end

  # Main traceline, if requested:
  if (!isscalar(lrgb))
    plot(Xdat,Y0+Ydat,"linewidth",lthk,"color",lrgb);
  end

end


######
##
function hln = VPline (xyxy,    # [x0 y0 x1 y1 ...] w.r.t. X0,Y0 baseline
                       ivp,     # VP index to anchor to
                       lspec)   # Line properties {prop, val, ...}
  global VP;
  X0 = VP.vpwidth(2);           # We assume reference is "left",
  Y0 = VP.lstY0(ivp);           # and "baseline".
  XX = X0 + xyxy(1:2:end);
  YY = Y0 + xyxy([2:2:end]);
  hln = line(XX, YY);
  set(hln, lspec{});  

end


######
##
function htxt = VPannotate (msg,        # Annotation text
                            ivp,        # VP index to anchor to
                            dx, dy,     # Offset from default position
                            side,       # Choose "left" or "right"
                            pos,        # Choose position option (see below)
                            tsize       # Font size *or* {prop, val, ...}
                           )
  global VP;

  if (isscalar(tsize))
    propval = {"fontsize", tsize};
  else
    propval = tsize;
  end

  X = VP.vpwidth(2);
  if (strcmp(tolower(side),"right"))
    X = VP.vpwidth(3);
  end

  Y = VP.lstY0(ivp);    # Default pos="base"
  valign = "bottom";
  if (strncmp(tolower(pos),"top",3))            # pos = "top"
     Y = VP.lstY0(ivp) + VP.lstDY(ivp);
     valign = "top";
  elseif (strncmp(tolower(pos),"above",5))      # pos = "above"
     Y = VP.lstY0(ivp) + VP.lstDY(ivp) + VP.vports{ivp}(2);
     valign = "top";
  elseif (strncmp(tolower(pos),"under",5))      # pos = "under"
     Y = VP.lstY0(ivp);
     valign = "top";
  elseif (strncmp(tolower(pos),"below",5))      # pos = "below"
     Y = VP.lstY0(ivp) - VP.vports{ivp}(3);
     valign = "bottom";
  elseif (strncmp(tolower(pos),"bottom",3))     # pos = "bot" (same as "base")
     Y = VP.lstY0(ivp);
     valign = "bottom";
  elseif (strncmp(tolower(pos),"base",4))       # pos = "base"
     # Just keep defaults here
  else
     printf("VPannotate: Unknown [pos]; assuming 'base'\n");
  end
  Y += dy;
  X += dx;

  htxt = text(X, Y, msg, "horizontalalignment", side,
              "verticalalignment", valign);
  set(htxt,propval{});
         
end

######
## Loads trace files into TR array
function GetTraces (tracefiles)

  if (!iscell(tracefiles))
    tracefiles = {tracefiles};
  end
  if (!iscellstr(tracefiles))
    error("Non-string argument for tracefiles.");
  end

  for i = 1:length(tracefiles)
    SEIS = load(tracefiles{i});
    if (isfield(SEIS,"STRACE"))        # Single scalar trace file...
      GetTraceFromSTRACEStruct(SEIS.STRACE);
    else                               # Or "Seismometer" file...      
      if (isfield(SEIS,"SEIS"))
        SEIS=SEIS.SEIS;
      end # Else assume old format not embedded in struct
      GetTracesFromSEISStruct(SEIS);
    end
  end

end

function GetTracesFromSEISStruct(SEIS)

      if (!isfield(SEIS,"AxesDesc"))
        SEIS.AxesDesc = "XYZ";
      end

      # Get unscaled signal maximum: (For plot range)
      ChFamMax_XYZ = max(max(SEIS.TraceXYZ));
      ChFamMax_PS  = max(max(SEIS.TracePS));
      ChFamMax_CPS = max(max(SEIS.CountPS));
      XRange = SEIS.TimeWindow;
      Tgpow = 2.0;  # Everything in trace file is proportional to
                    # energy (Tgpow=2.0). At least for current SIES
                    # format specification.

      ChanMax = max(SEIS.TraceXYZ(:,1));    ## Trace X
      ChLabel = SEIS.AxesDesc(1);
      Data = SEIS.TraceXYZ(:,1);
      AddTrace(Data, XRange, ChanMax, ChFamMax_XYZ, Tgpow, ChLabel, "");

      ChanMax = max(SEIS.TraceXYZ(:,2));    ## Trace Y
      ChLabel = SEIS.AxesDesc(2);
      Data = SEIS.TraceXYZ(:,2);
      AddTrace(Data, XRange, ChanMax, ChFamMax_XYZ, Tgpow, ChLabel, "");

      ChanMax = max(SEIS.TraceXYZ(:,3));    ## Trace Z
      ChLabel = SEIS.AxesDesc(3);
      Data = SEIS.TraceXYZ(:,3);
      AddTrace(Data, XRange, ChanMax, ChFamMax_XYZ, Tgpow, ChLabel, "");

      ChanMax = max(SEIS.TracePS(:,1));     ## Trace P
      Data = SEIS.TracePS(:,1);
      AddTrace(Data, XRange, ChanMax, ChFamMax_PS, Tgpow, "P", "");

      ChanMax = max(SEIS.TracePS(:,2));     ## Trace S
      Data = SEIS.TracePS(:,2);
      AddTrace(Data, XRange, ChanMax, ChFamMax_PS, Tgpow, "S", "");

      ChanMax = max(SEIS.CountPS(:,1));     ## Count P
      Data = SEIS.CountPS(:,1);
      AddTrace(Data, XRange, ChanMax, ChFamMax_CPS, Tgpow, "p", "");

      ChanMax = max(SEIS.CountPS(:,2));     ## Count S
      Data = SEIS.CountPS(:,2);
      AddTrace(Data, XRange, ChanMax, ChFamMax_CPS, Tgpow, "s", "");

end

function GetTraceFromSTRACEStruct(TRACE)

      XRange = TRACE.TimeWindow - TRACE.OTime;  # W.r.t event origin time

      if (isfield(TRACE,"Tgpow"))
         Tgpow = TRACE.Tgpow;
      else                      # Else assume trace contains
        Tgpow = 1.0;            # amplitude data
      end

      ChanMax = max(TRACE.Trace(:,1));
      ChLabel = "";
      Data = TRACE.Trace(:,1);
      AddTrace(Data, XRange, ChanMax, ChanMax, Tgpow, ChLabel, "");
      AddTrace(-Data, XRange, ChanMax, ChanMax, Tgpow, "N", "");

end

function AddTrace(Data, XRange, ChanMax, ChFamMax, Tgpow, ChLabel, ChDesc)
  global TR;
  TR{end+1}.Data = Data;        # [y1 y2 y3 ...]
  TR{end}.XRange = XRange;      # [xmin xmax]
  TR{end}.ChanMax = ChanMax;    # Channel max value
  TR{end}.ChFamMax = ChFamMax;  # Max value among family of channels
  TR{end}.Tgpow = Tgpow;        # Trace gpow (1.0: trace is amp, 2.0: energy)
  TR{end}.ChLabel = ChLabel;    # Label for tick mark ('R' for Radial, eg)
  TR{end}.ChDesc = ChDesc;      # Brief channel description

  if true
    printf("%s %2i: '%1s'. %s %8.4e (%8.4e), %s: %8.4e, %s %4.2f.\n",
           "Adding channel", length(TR), ChLabel, 
           "ChanMax:", ChanMax, sqrt(ChanMax), 
           "ChFamMax:", ChFamMax,
           "Tgpow:", Tgpow);
  end
end

######
## If we receive a vector, assume explicit color vector and return it.
## If we receive a scalar >= 1, assume index into colormap and return color.
## If we receive a scalar < 1, return 0 as a non-color signal
function color = ColorMap(idx)

  if (!isscalar(idx))
     color = idx;
     return
  end

  map = [ [0.85 0.00 0.00];     # 1 color_x   (red)
          [0.00 0.60 0.00];     # 2 color_y   (green)
          [0.00 0.10 0.75];     # 3 color_z   (blue)
          [0.90 0.00 0.60];     # 4 color_p   (pink)
          [0.00 0.70 0.90];     # 5 color_s   (cyan)
          [1.00 0.70 0.90];     # 6 color_pf  (light pink)
          [0.70 0.90 1.00];     # 7 color_sf  (light cyan)
          [0    0    0   ];     # 8 black
          [1.00 0.60 0.20];     # 9 pick_1    (orange)
          [1.00 0.60 0.20];     # 10 pick_2   (orange)
          [1.00 0.60 0.20];     # 11 pick_3   (orange)
          [1.00 0.60 0.20];     # 12 pick_4   (orange)          
        ];

  color = 0;
  idx = floor(idx);
  register = floor((idx-1)/length(map));
  if (idx>0);
    idx=1+mod(idx-1,length(map)); 
    color = map(idx,:);
  end
  if (register>0);  # "Secret" colors; (subject to change w/o notice)
    tgtcolor = [1 1 1];
    cw=1;
    if (register>10);
       tgtcolor = [0 0 0];
       register -= 10;
       cw=2;  # Extra weight to chosen color; shades are harsh.
    end
    color = ((cw*color + register*tgtcolor) # (each "register" gets
             / (register+cw));              #  us closer to white,
  end                                       #  or black.)
end
