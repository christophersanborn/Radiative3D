# normcurve_compare.m
#
# Creates a plot of two NORMCURVE series representing integrated
# energy vs. time for two different runs, alongside a third NORMCURVE
# series representing a reference curve, generally a best-fit curve to
# the baseline curve.
#
# Note: Plot visibility is set to "off" by default to speed up
# scripted runs.  If making plot interactiviely, may need to
# set(gcf(),"visible", "on") in order to see it.
#

function normcurve_compare(
                NC1,              # Either a NORMCURVE struct or a cell array
                                  # of NORMCURVE structs
                NC2,              # Either a NORMCURVE struct or [], if NC1 is
                                  # a cell array.
                REF,              # NORMCURVE for Reference curve
                clipindex=20,     # Ignore indices before this when computing
                                  # view window y-max.
                ROI={},           # Region of Interest struct or array of
                                  # structs
                ViewY = [0 0],       # Override Y bounds of view window
                                     # (if non-zero)
                PaperWH = [6.0 4.5], # Paper width, height

                RangePower = 0    # Scale energy by r^q where q is RangePower
                                  # Can be used to remove geometric spreading.
           )

  # Default Color sequence:
  Colors{1} = [0.1 0.1 0.7];  # Blue
  Colors{2} = [0.8 0.4 0.1];  # Orange
  Colors{3} = [0.1 0.7 0.1];  # Green
  Colors{4} = [0.7 0.1 0.7];  # Purple
  Colors{5} = [0.7 0.1 0.1];  # Red (Reference color)

  # Domain:
  X = linspace(REF.RangeWindow(1), REF.RangeWindow(2),
               length(REF.SummedEnergy));

  # Series:
  R1 = REF.SummedEnergy;
  R1txt="";
  if (isfield(REF,"Label")); R1txt = REF.Label; end
  if (!isfield(REF,"Visible")); REF.Visible = true; end
  # Unpack args into a cell array of Y-series.
  if (iscell(NC1))
    for i=1:length(NC1)
      Y{i} = NC1{i}.SummedEnergy;
      Ytxt{i} = "";
      Ycolor{i} = Colors{mod(i-1,length(Colors))+1};
      Ylinestyle{i} = "-";
      if (isfield(NC1{i},"Label")); Ytxt{i} = NC1{i}.Label; end
      if (isfield(NC1{i},"ColorIndex")); Ycolor{i} = Colors{NC1{i}.ColorIndex}; end
      if (isfield(NC1{i},"LineStyle")); Ylinestyle{i} = NC1{i}.LineStyle; end
    end
  else
    Y{1} = NC1.SummedEnergy;
    Y{2} = NC2.SummedEnergy;
    Ytxt{1} = Ytxt{2} = "";
    Ylinestyle{1} = Ylinestyle{2} = "-";
    Ycolor{1} = Colors{1};
    Ycolor{2} = Colors{2};
    if (isfield(NC1,"Label")); Ytxt{1} = NC1.Label; end
    if (isfield(NC1,"ColorIndex")); Ycolor{1} = Colors{NC1.ColorIndex}; end
    if (isfield(NC1,"LineStyle")); Ylinestyle{1} = NC1.LineStyle; end
    if (isfield(NC2,"Label")); Ytxt{2} = NC2.Label; end
    if (isfield(NC2,"ColorIndex")); Ycolor{2} = Colors{NC2.ColorIndex}; end
    if (isfield(NC2,"LineStyle")); Ylinestyle{2} = NC2.LineStyle; end
  end

  # Normalize on end-value of Y{1}, (Assumes Y{1} represents the "baseline"
  # condition, so all others will compare with that), and scale enery by
  # r^q, e.g. to remove effects of geometric spreading. (We pick unity as
  # the end point so that Y{1}(end) will still be 1.0.)
  E_baseline = Y{1}(end) * (X(end) ./ X').^RangePower;
  R1 = R1 ./ E_baseline;
  for i=1:length(Y)
      Y{i} = Y{i} ./ E_baseline;
      smallest(i) = min(Y{i});  # Used in view-window determination
  end

  # Determine view window:
  clipindex = 20;
  vymax = max(Y{1}(clipindex:end));
  vymin = min(smallest);
  viewwindow = [X(1), X(end), vymin/2, vymax];
  printf("Detected view window [xmin xmax ymin ymax] = [%g %g %g %g]\n",
         viewwindow);
  if (sum(ViewY.^2)>0)  # Then override Y bounds
    viewwindow = [X(1), X(end), ViewY(1), ViewY(2)];
    printf("Overriding view wndw [xmin xmax ymin ymax] = [%g %g %g %g] %s\n",
           viewwindow, "at user's request.");
  end
  viewwindow += [1 0 0 0]; # (Prevent 0 from appearing on axis.)
                           # ((Mild aesthetic preference...))

  ##
  ## ...And PLOT!
  ##

  txtsize = 9.0;                # Axes ticktext and labels
  titlesize = 9.0;              # Title text slightly bigger
  legtxtsize = 10.0;            # And legend text slightly smaller
  serieslinewidth = 8;          # For the E-curve series
  thinlinewidth = 2.0;          # For the reference curve

  figinit(PaperWH(1), PaperWH(2),   # Clear and initialize 
      "paperunits", "inches",       # a 4.0" x 3.0" figure
      "visible", "off");            #
  axes("nextplot", "add",           # Create and initialize an axes object
       "fontweight", "bold", "fontsize", txtsize, "linewidth", 1.0);

  hNoLine = semilogy(1,1,"color", "white", "linestyle", ":"); # Create a gap
  legend(hNoLine," ");  # Legend Null-Entry; Keeps other entries out of corner

  printf("  Ref curve: E(x)=c*x^q; c=%g, q=%g\n",...
                    REF.PLCQ_Summed(1),REF.PLCQ_Summed(2));
  if (REF.Visible)  # Reference curve
    hR1 = semilogy(X, R1, "color", Colors{end},
                      "linewidth", thinlinewidth, "linestyle", "--");
    if (length(R1txt)>0)
      legend(hR1, R1txt);
    end
  end

  for i=1:length(Y)             # For each e-curve series:
    LWidth = serieslinewidth;   # default linewidth, we might adjust per-series
    if (strcmp(Ylinestyle{i},":")); LWidth *= 1.5; end;
    hY{i} = semilogy(X, Y{i}, "color", Ycolor{i},
                     "linestyle", Ylinestyle{i}, "linewidth", LWidth);
    if (length(Ytxt{i})>0)
      hleg = legend(hY{i}, Ytxt{i});
    end
    printf("  Series #%i reduction (dB): %14g  \"%s\"\n", 
           i, 10*log10(Y{i}(end)), Ytxt{i});
  end

  semilogy(X, Y{1}, "color", Ycolor{1},   # Plot baseline again so it layers on
           "linewidth", serieslinewidth,  # top (without messing up legend order)
           "linestyle", Ylinestyle{1});   #

  axis(viewwindow);
  annotate_regions(ROI);

  title ("Energy Comparison, Whole Curve", "fontsize", titlesize);
  xlabel("Range (km)");
  ylabel(sprintf("Sum(Energy*dt) (Relative)"));
  set(gca(), "yaxislocation","right",
             "yminortick", "on", "yminorgrid", "off"#,
            #"ticklength", [0.015 0.025]  # 50% bigger tick size
     );

  hlgnd = findobj(gcf(),"type","axes","Tag","legend");  # Set legend text size
  set(hlgnd, "fontweight", "normal", "FontSize", legtxtsize, "FontName", "Mono");

end

######
# Helper Functions:
#

