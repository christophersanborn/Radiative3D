## fehlerR1R2plot.m
##
##
function fehlerR1R2plot(
                R1R2data
           )

  B0 = R1R2data(:,1);
  QS = R1R2data(:,2);
  QI = R1R2data(:,3);
  R1 = R1R2data(:,4);
  R2 = R1R2data(:,5);

  [uniqueQS, dummy, uniqueQSidx] = unique(QS);
  [uniqueQI, dummy, uniqueQIidx] = unique(QI);

  figinit(5.0, 4.5,#3.75,                # Clear and initialize
       "paperunits", "inches",      # a 5.0" x 3.75" figure
       "visible", "off");           #
  axes("fontsize", 9,               # Create and initialize axes obj
       "linewidth", 0.5, "nextplot", "add");

  QSColors = {[.85 0 0], [0 .75 0], [0 0 .85]};
  QIColors = {[1 1 1]*.65, [1 1 1]*.35, [1 1 1]*0};

  # Segregate on QS, Sort on QI
  Colors = QSColors;
  SegIdxSet = uniqueQSidx;
  SortSet = QI;
  for (i=1:length(uniqueQS))
    uidx = find(SegIdxSet==i);
    [dummy, sortidx] = sort(SortSet(uidx));
    suidx = uidx(sortidx);
    R1set = R1(suidx);
    R2set = R2(suidx);
    h = plot(R1set,R2set, "color", Colors{mod(i-1,length(Colors))+1},
                          "linewidth", 5
        );
    ltext = sprintf("Q_{Scat} = %7.2f", uniqueQS(i));
    legend([h],ltext, "location", "northwest");
    legend("boxoff");
  end

  # Segregate on QI, Sort on QS
  Colors = QIColors;
  SegIdxSet = uniqueQIidx;
  SortSet = QS;
  for (i=1:length(uniqueQI))
    uidx = find(SegIdxSet==i);
    [dummy, sortidx] = sort(SortSet(uidx));
    suidx = uidx(sortidx);
    R1set = R1(suidx);
    R2set = R2(suidx);
    h = plot(R1set,R2set, "color", Colors{mod(i-1,length(Colors))+1},
                          "linewidth", 2.0
        );
    ltext = sprintf("Q_{Int} =\t%4g", uniqueQI(i));
    hl = legend([h],ltext, "location", "northwest");
    set(hl, "fontsize", 7, "fontname", "Mono");
    legend("boxoff");
  end

  get(legend())

  # Data Labels
  for (i=1:length(R1))
    text(R1(i), R2(i), sprintf("%d",i),
         "horizontalalignment", "right",
         "verticalalignment", "bottom",
         "fontname", "Courier", "fontsize", 8);
    text(R1(i), R2(i), sprintf("%.2f",B0(i)),
         "horizontalalignment", "left",
         "verticalalignment", "top",
         "fontname", "Times", "fontsize", 6);
  end

  # Labels and Titles
  labelfontname = "Helvetica";
  labelfontsize = 7;
  xlabel("R1 = log10[E_{early} @ r = 150 km / E_{late} @ r = 150 km]",
        "fontweight", "bold", "fontsize", labelfontsize, "fontname", labelfontname);
  ylabel("R2 = log10[r^2E_{early} @ r = 50 km / r^2E_{early} @ r = 150 km]",
        "fontweight", "bold", "fontsize", labelfontsize, "fontname", labelfontname);
  title("R_2 vs. R_1 in Preliminary Model", "fontweight", "bold",
                    "fontsize", labelfontsize+3, "fontname", labelfontname);

end
