## arrayimage.m
##
## Plots a traveltime curve from a seismic array using image() plot.
##
##
## Printing to file:
##
##   print("fig.png", "-r160");  # for  800x600 png file
##   print("fig.png", "-r240");  # for 1200x900 png file
##
##   Do not use the "-S" option to print - it won't do what you want.
##
## Common Normalization:
##
##  The return value, NS, is a struct containing energy curves that
##  can be processed and passed back on a subsequent run to control
##  image normalization.
##
## Note: Plot visibility is set to "off" by default to speed up
## scripted runs.  If making plot interactiviely, may need to
## set(gcf(),"visible", "on") after this function in order to see it.
##
function NS = arrayimage (           # NS is "Norm Struct"
                  ARRAY,             # An ARRAY struct (see array.m)
                  AXES = [1 1 1],    # Which axes to include
                  gamma = 1.0,       # gamma factor for mag scaling
                  timewindow = [],   # time window for viewing
                  norm = 0.3,        # Norm ratio (legacy) or NORMCURVE struct:
                                     #  1.0: peak-based normalization
                                     #  0.0: area-based normalization
                  ColorLim = 0       # Upper limit of color axis (0 for auto).
                )

  # Look for MParams file:
  if (exist("out_mparams.octv","file")==2)
    MPAR=load("out_mparams.octv");
    NPhonCast = MPAR.NumPhonons;
  else
    NPhonCast = -1;
  end

  BB = arraymatrix(ARRAY, [1 1 1]);  # Get matrix form of trace

  dt = (ARRAY.TimeWindow(2)-ARRAY.TimeWindow(1)) / ARRAY.NumBins;
  EE = sum(BB,2) * dt;               # Get summed energy by seismometer.
                                     # (Will be plotted as an annotation)

  ## Assemble Norm Struct for output:
  NS.TimeWindow = ARRAY.TimeWindow;
  NS.RangeWindow = [ARRAY.Distances(1), ARRAY.Distances(end)];
  NS.SummedEnergy = EE;
  NS.PeakEnergy = max(BB,[],2);

  if (isstruct(norm)) ## Then normalize via NORMCURVE struct:

    NormCurve = norm.SummedEnergy / (norm.TimeWindow(2)-norm.TimeWindow(1));
    BB = BB ./ NormCurve;   # Elements are now fraction of time-
                            # averaged curve value
    BB = BB .^ (1/gamma);   # Gamma-scaling

  else  ## Else legacy normalization

    #  Apply gamma scaling  (Helps when the tall peaks make the small
    #                        peaks hard to see)
    BB = BB .^ (1/gamma);   # A gamma factor of 2.0 gives us square-root
                            # scaling

    #  Individually normalize each row:
    #                       (Helps prevent more distant traces from
    #                       washing out due to their lower overall
    #                       amplitude compared with the near traces.)

    BB_area = BB ./ sum(BB,2);    ## (Matrix ./ works with Octave >= 3.6.0,
                                  ##  otherwise produces an error)
                                  ## Normalizes each row based on area 
                                  ## under curve.

    BB_peak = BB ./ max(BB,[],2); ## Normalizes each row based on 
                                  ## peak value.

    BB = (1-norm)*BB_area + norm*BB_peak;   ## Splitting the difference seems
                                            ## to give nicest overall effect.
  end

  ## Physical dimensions

  SamplesPerSeis = length(ARRAY.TraceXYZ) / ARRAY.NumSeismometers;
  X0 = ARRAY.Distances(1);
  X1 = ARRAY.Distances(end);
  X = linspace(X0, X1, ARRAY.NumSeismometers);
  T0 = ARRAY.TimeWindow(1);
  T1 = ARRAY.TimeWindow(2);
  T = linspace(T0, T1, SamplesPerSeis);

  if numel(timewindow) != 2
    timewindow = [T0 T1]; # default viewing window = data window
    # else time window was given in args
  end

  ## Plot
  #
  #  plots with time on Y axis, dist on X.  This seems to be the
  #  standard for TT graphs.  Still, may want to reverse this as it
  #  seems somewhat foreign to me.
  #

  linethin   = 0.5;
  linethick1 = 3.75;
  fonttitle  = 9.5;
  fontxylabel = 9.0;
  fontaxes = 9.0;

  figinit(5.0, 3.75,                # Clear and initialize 
       "paperunits", "inches",      # a 5.0" x 3.75" figure
       "visible", "off");           #
  axes("fontsize", fontaxes,        # Create and initialize axes obj
       "linewidth", linethin, "nextplot", "add");

  imagesc(X, T, BB');
  axis("xy");
  axis([X0 X1 timewindow(1) timewindow(2)]);
  if (ColorLim>0); caxis([0 ColorLim]); end
  map = hot(64);
  map = map(end:-1:1,:);  # light to dark instead of dark to light
  colormap(map);
  colorbar("linewidth", linethin, "fontsize", fontxylabel, "eastoutside");
  xlabel("Range (km)", "fontsize", fontxylabel);
  ylabel("Time (s)", "fontsize", fontxylabel);
  azi_deg = ARRAY.Azimuths(end);
  title(sprintf("Travel-time curve at %g-degrees azimuth", azi_deg),
        "fontsize", fonttitle, "fontweight", "bold");
 
  # Overlays:
  RefCurve = [];
  RefCaption = "";
  OLTextSize = 8.5;   # Overlay text size
  OLY0 = 0.80;      # Base and height of refcurve overlay 
  OLdY = 0.15;      #            *
  CaptLR = "";      # Low-right caption
  OLCscale = 0.5;   # Overlay Curve scale (0.5 sqr rt)
  CaptUR = "Scale: Amplitude";
  if (isstruct(norm)) # Extra goodies if norm is a struct:
    RefCurve = norm.SummedEnergy;
    if (isfield(norm,"PLCQ_Summed"))
      RefCaption = sprintf(
                    "\n\nRef Curve: E(x)=c*x**q;\nc = %8.2e;\nq = %0.2f;",
                    norm.PLCQ_Summed(1), norm.PLCQ_Summed(2));
      OLTextSize = 7; # (smaller text for extra caption content)
      OLY0 = 0.825;   # (smaller ref-curve vertical size)
      OLdY = 0.125;   #
    end
    OLCscale = 0.25;
    CaptUR = "Curve Scale: 4th-root";
  else
    CaptLR = sprintf("NR: %0.2f", norm);
  end
  CutIn = max(ceil(length(EE)/10),2);
  EEsubset = EE(CutIn:end); # Blows up near zero, so excise first elements
  Xsubset = X(CutIn:end);   #
  [PeakE,iPE] = max(EEsubset);  # Max value of E (after cut-on)
  TermE = EEsubset(end);        # Terminal value of E
  explfrac = getexplfrac(ARRAY.EventMT);
  caption = sprintf("Sum(E*dt);  Max: %8.2e at %1.0f km;%s%s%s%s%s%s", 
                    PeakE, Xsubset(iPE),
                    sprintf("\n           Term: %8.2e at %1.0f km;",
                            EEsubset(end), Xsubset(end)),
                    sprintf("\n\nPhonons Cast: %8.2e", NPhonCast),
                    sprintf("\nEvent Isotropy: %0.0f%%", explfrac*100),
                    sprintf("\nFrequency: %0.2f Hz", ARRAY.Frequency),
                    RefCaption);

  overlayplot(OLY0, OLdY, X(CutIn:end),
              EE(CutIn:end).^OLCscale, RefCurve(CutIn:end).^OLCscale,
              {caption, CaptUR, CaptLR}, OLTextSize);

end #END FUNCTION

######
## FUNC: overlayplot()
##   Draws a nice set of overlays with a trace of a data series
##   provided, and up to three caption positions.
##
function overlayplot (Y0,       # Y-origin, as a fraction of view window
                                #  (ie, 0.5 is middle of plot)
                      dY,       # Peak height, as fraction of view window.
                      X, Y,     # X and Y vectors. Y will be scaled to fit
                                #  within Y0+/-dY.
                      R,        # A Y-series that serves as a reference
                                # curve (pass empty matrix [] if none.)
                      caption,  # Caption text or caption array
                      tsize=7,  # Text size
                      lwidth=2  # Base line width
                     )

  # Up to three captions: (main, upper-right, lower-right)
  caption3 = caption2 = "";
  if (iscell(caption))
    captions = {caption{},"",""}; # (ensure length>=3)
    caption1 = captions{1};
    caption2 = captions{2};
    caption3 = captions{3};
  else
    caption1 = caption;
  end

  # View window dimensions:
  axlims = axis();
  viewheight = axlims(4)-axlims(3);
  viewwidth  = axlims(2)-axlims(1);
  plot_y0 = axlims(3) + Y0*viewheight;
  plot_dy = dY*viewheight;
  rightmarg = axlims(1)+0.96*viewwidth;
  topmarg = axlims(3)+0.96*viewheight;
  botmarg =  axlims(3)+0.05*viewheight;

  # Scale curves to view window:
  [PeakY, iPY] = max(abs(Y));
  Yb = plot_y0 + 0*Y;
  YY = plot_y0 + plot_dy*(Y/PeakY);
  RR = plot_y0 + plot_dy*(R/PeakY);

  # Plot curves:
  lwidth_th = lwidth/1.5;
  tag="OLCurve";
  if (length(RR)>0) # Reference curve, if present
     plot(X, RR, "color", "black", "linestyle", ":",
          "linewidth", lwidth_th, "tag", tag);
  end
  plot(X, Yb, "color", "black", "linewidth", lwidth, "tag",tag); # Flatline
  plot(X, YY, "color", "black", "linewidth", lwidth, "tag",tag); # Scaled trace
  plot([X([iPY,iPY])],[plot_y0, YY(iPY)], "color", "black", # Peak mark
       "linestyle", ":", "linewidth", lwidth_th, "tag", tag);

  # Captions:
  text(X(1), plot_y0, caption1, "fontsize", tsize,  # Main caption
       "horizontalalignment", "left",
       "verticalalignment", "top", "tag", "OLCapt1");
  text(rightmarg, topmarg, caption2, "fontsize", tsize,  # Upper-right caption
       "horizontalalignment", "right",
       "verticalalignment", "top", "tag", "OLCapt2");
  text(rightmarg, botmarg, caption3, "fontsize", tsize,  # Low-right caption
       "horizontalalignment", "right",
       "verticalalignment", "bottom", "tag", "OLCapt3");
      # NOTE: The captions are "tagged", to facilitate, e.g., knocking out
      # a particular caption for publication, with something like:
      #     set(findobj("tag", "OLCapt3"), "color", [1 1 1]);

end

function expl = getexplfrac(MT)

  tr = trace(MT);

  explpart = (tr^2) / 3;   ## Mag^2 of explosive part
  total = sum(sum(MT.*MT)); ## Mag^2 of whole tensor

  expl = explpart/total;
  if tr<0; expl*=-1; end

end
