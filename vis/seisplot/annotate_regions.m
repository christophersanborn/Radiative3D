# annotate_regions.m
#
function annotate_regions(ROI)

  if (!iscell(ROI))     # If not cell array...
    RA{1} = ROI;        # then make it a cell array
  else                  #
    RA = ROI;
  end

  VW = axis();  # Get the view window

  logscale = strcmp("log", get(gca(),"yscale"));
  if (logscale)
    bot = log(VW(3));  # Used to compute paper positions
    top = log(VW(4));  #
    height = top-bot;
  else
    bot = VW(3);
    top = VW(4);
    height = top-bot;
  end

  for i=1:length(RA)
    Xwin = RA{i}.RegionWindow;
    Xmid = 0.5*(sum(Xwin));
    Yleft(1) = bot + RA{i}.RegionVSpan(1)*height;
    Yleft(2) = bot + RA{i}.RegionVSpan(2)*height;
    Ytxt = bot + RA{i}.LabelVPos*height;
    if (logscale)
      Yleft = exp(Yleft);
      Ytxt = exp(Ytxt);
    end
    Yright = Yleft;
    plot(Xwin(1)*[1 1], Yleft, "color", "black",
                               "linewidth", RA{i}.RegionLineWidth,
                               "linestyle", ":");
    plot(Xwin(2)*[1 1], Yright, "color", "black",
                                "linewidth", RA{i}.RegionLineWidth,
                                "linestyle", ":");
    text(Xmid,Ytxt,RA{i}.Label,
           "horizontalalignment", "center",
           "verticalalignment", "middle", "fontsize", RA{i}.LabelFontSize);

  end

end
