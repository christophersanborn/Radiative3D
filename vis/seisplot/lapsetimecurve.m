## arrayimage.m
##
## Plot a lapse-time curve on existing plot axes
##
function NS = lapsetimecurve (            # NS is "Norm Struct"
                  ARRAY,                  # An ARRAY struct (see array.m)
                  QScatInt,               # Scattering and Intrinsic Q ([q1 q2])
                  AXES = [0 0 1],         # Which axes to include
                  geospread = 2,          # Geometric spreading exponent
                  PhaseEdge = [3.6, 0],   # v, t0 for phase
                  LWindow1 = [5 20],      # Lapse window 1
                  LWindow2 = [45 115]     # Lapse window 2
                )

  # Look for MParams file:
  if (exist("out_mparams.octv","file")==2)
    MPAR=load("out_mparams.octv");
    NPhonCast = MPAR.NumPhonons;
  else
    NPhonCast = -1;
  end

  QScat = QScatInt(1);
  QInt = QScatInt(2);
  Albedo = QScat^-1 / (QScat^-1 + QInt^-1);

  [dum, idx8]   = min(abs(ARRAY.Distances - 8));    # index closest to 8km
  [dum, idx50]  = min(abs(ARRAY.Distances - 50));   # index closest to 50km
  [dum, idx150] = min(abs(ARRAY.Distances - 150));  # index closest to 150km
  printf("Distances %f, %f\n", ARRAY.Distances(idx50), ARRAY.Distances(idx150));

  BB = arraymatrix(ARRAY, AXES);     # Get matrix form of trace

  dt = (ARRAY.TimeWindow(2)-ARRAY.TimeWindow(1)) / ARRAY.NumBins;
  EE = sum(BB,2) * dt;               # Get summed energy by seismometer.
                                     # (Will be plotted as an annotation)

  ## Lapse-Window Integrals: Window 1:
  EE1 = REE1 = zeros(ARRAY.NumSeismometers,1);  # Sums E*dt and R^2E*dt
  Times1 = zeros(ARRAY.NumSeismometers,2);
  for iseis = 1:ARRAY.NumSeismometers
    t_edge = PhaseEdge(2) + ARRAY.Distances(iseis)/PhaseEdge(1);  # t0+r/v
    t_offset = LWindow1(1);
    t_duration = LWindow1(2)-LWindow1(1);
    t_begin = t_edge + t_offset;
    iwinbegin = max(1, ceil(t_begin/dt));
    iwinend = iwinbegin + round(t_duration/dt) - 1;
                #              !-----!-----!-----!-----!-----!
                #  |     |     |     |     |     |     |     |     |     |     |
                #     1     2     3     4     5     6     7     8     9     0     1
                #                 B                       E
    Times1(iseis,:) = [(iwinbegin-1) iwinend]*dt;
    tempEE = BB(iseis,iwinbegin:iwinend);
    EE1(iseis) = sum(tempEE)*dt;
    REE1(iseis) = sum(tempEE)*dt * (ARRAY.Distances(iseis))^geospread;
  end
  ## Lapse-Window Integrals: Window 2:
  EE2 = REE2 = zeros(ARRAY.NumSeismometers,1);  # Sums E*dt and R^2E*dt
  Times2 = zeros(ARRAY.NumSeismometers,2);
  for iseis = 1:ARRAY.NumSeismometers
    t_edge = PhaseEdge(2) + ARRAY.Distances(iseis)/PhaseEdge(1);  # t0+r/v
    t_offset = LWindow2(1);
    t_duration = LWindow2(2)-LWindow2(1);
    t_begin = t_edge + t_offset;
    iwinbegin = ceil(t_begin/dt);
    iwinend = iwinbegin + round(t_duration/dt) - 1;
                #              !-----!-----!-----!-----!-----!
                #  |     |     |     |     |     |     |     |     |     |     |
                #     1     2     3     4     5     6     7     8     9     0     1
                #                 B                       E
    Times2(iseis,:) = [(iwinbegin-1) iwinend]*dt;
    tempEE = BB(iseis,iwinbegin:iwinend);
    EE2(iseis) = sum(tempEE)*dt;
    REE2(iseis) = sum(tempEE)*dt * (ARRAY.Distances(iseis))^geospread;
  end

  if (false) # print time windows
    iseis=[1,8:8:length(ARRAY.Distances)];
    [ARRAY.Distances(iseis) Times1(iseis,:) Times2(iseis,:)]
  end

  global linestyle;

  hold on;

  refval = REE1(idx8);  # Reference value against which curves are plotted

  h1 = semilogy(ARRAY.Distances,REE1/refval,"color", [0 0 .6],
           "linewidth", linestyle.width, "linestyle", linestyle.style);
  h2 = semilogy(ARRAY.Distances,REE2/refval,"color", [0 .6 0],
           "linewidth", linestyle.width, "linestyle", linestyle.style);

  FehlerR1 = log10(EE1(idx150)/EE2(idx150));
  FehlerR2 = log10(REE1(idx50)/REE1(idx150));

  ##LegendTxt = sprintf("B_0, Q_{scat}, Q_{int} = %6.2f, %6.0f, %6.0f",
  LegendTxt = sprintf("%4.2f  %7.0f %7.0f  %7.4f  %7.4f",
                      Albedo, QScat, QInt, FehlerR1, FehlerR2);
  legend([h1],LegendTxt);
  legend("boxoff");

end


######
## Perhaps there's an easier way to do this, but... appends legend
## item, preserving forgoing items, and specifies style.
##
## Ack, basically none of this is needed, lol.  If you call legend
## with an array of line handles (eg. from h = plot(...);
## legend(h,"leg text");) then it already appends to the legend, and
## does exactly what I want.  (I can even fake it by doing
## hfake=plot(NaN,NaN,propvals...) if I want style to differ from
## actual displayed traces.)  This solve the problem of (1) appending,
## and (2) preventing the persistent and pernicious behavior of taking
## the props from lines that are not intended for the legend.
##
## So... don't use this anymore, lol. (And it didn't work anyway.)
##
function hleg = LegendAppend(txt, color, thick, style="-")

  nitems = 0;

  if (legend()) hleg = legend(); else hleg=0; end
  if (hleg) # Get existing items
    hchildren = get(hleg, "children");
    nitems = length(hchildren)/2;
    hlines = hchildren(1:nitems);
    htxts = hchildren((nitems+1):end);
    for i = 1:nitems
      artxt{end+1} = get(htxts(i),"string");
      arcolor{end+1} = get(hlines(i),"color");
      arthick{end+1} = get(hlines(i),"linewidth");
      arstyle{end+1} = get(hlines(i),"linestyle");
    end
  end

  artxt{end+1} = txt;
  arcolor{end+1} = color;
  arthick{end+1} = thick;
  arstyle{end+1} = style;
  nitems += 1;

  hleg = legend(artxt);
  hchildren = get(hleg, "children");
  hlines = hchildren(1:nitems);
  htxts = hchildren((nitems+1):end);
  for i = 1:nitems
    set(hlines(i), "color", arcolor{i});
    set(hlines(i), "linewidth", arthick{i});
    set(hlines(i), "linestyle", arstyle{i});
  end

end



### TODO discard below, mostly
function rando()

#########
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
  colorbar("eastoutside", "linewidth", linethin, "fontsize", fontxylabel);
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
              EE(CutIn:end).^OLCscale, RefCurve(16:end).^OLCscale, 
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
