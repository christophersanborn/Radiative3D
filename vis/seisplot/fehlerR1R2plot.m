## fehlerR1R2plot.m
##
##
function fehlerR1R2plot(
                R1R2data,
                PaperWH = [5.0 4.5],    # Paper width, height
                boolDataLabels = true,  # To number R1R2 points or not
                distandspread = [50 150 2], # Near, Far, and geo spread exponent
                freq = 0,
                veloc = 0,
                R1params = [0.1 30],    # [late-fraction, late-t-offset]
                FehlerSatoRatio = [],   # If non-empty/zero, fixed Qfehler multiple
                                        # for legend text.
                ScatLegLoc = "northwest",
                IntLegLoc = ScatLegLoc  # (Don't set they don't separate)
           )

  B0 = R1R2data(:,1);
  QS = R1R2data(:,2);
  QI = R1R2data(:,3);
  R1 = R1R2data(:,4);
  R2 = R1R2data(:,5);

  RangeNear = distandspread(1);
  RangeFar  = distandspread(2);
  GeoSpread = distandspread(3);
  R1alpha = R1params(1);
  LateWindowToffset = R1params(2);
  if (length(FehlerSatoRatio)!=1); FehlerSatoRatio=0; end

  [uniqueQS, dummy, uniqueQSidx] = unique(QS);
  [uniqueQI, dummy, uniqueQIidx] = unique(QI);

  figinit(PaperWH(1), PaperWH(2),   # Clear and initialize
       "paperunits", "inches",      # a 5.0" x 4.5" figure
       "visible", "off");           #
  axes("fontsize", 9,               # Create and initialize axes obj
       "linewidth", 0.5, "nextplot", "add");

  QSColors = {[.8 0 .6], [.85 0 0], [0.8 0.4 0], [.82 .77 0], [0 .65 0],...
              [0 0.45 0.65], [.1 .0 .85], [.65 0 .75]};
  QIColors = {[1 1 1]*.75, [1 1 1]*.65, [1 1 1]*.60, [1 1 1]*.55, [1 1 1]*.45,...
              [1 1 1]*.35, [1 1 1]*.25, [1 1 1]*.15, [1 1 1]*0,   [1 1 1]*0};

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
    ltext = sprintf("Q_{Sato} = %7.2f", uniqueQS(i));
    legend([h],ltext, "location", ScatLegLoc);
    legend("boxoff");
    # (Compute measured Scat Q's)
    if (freq != 0)
      R2mult = 2.302585 * veloc / (2*pi*freq*(RangeFar-RangeNear));
      QscatR2List = [];
      QscatR1List = [];
      table = [];
      for (j = suidx(1:end-2)')
        # Solve Qscat from known Qint and R2
        Qtot = 1/(R2(j) * R2mult);
        QscatR2 = 1/(Qtot^-1 - 1/QI(j));
        QscatR2List(end+1) = QscatR2;
        # Solve Qscat from known R1 and R1org
        R1org = log10(R1alpha);
        R1mult = 2*pi*freq*LateWindowToffset;
        QscatR1 = (2*pi*freq*RangeFar/veloc) / ...
                  log(1 + exp(R1mult/QI(j) - log(10)*(R1(j)+R1org)));
        QscatR1List(end+1) = QscatR1;

        qtext = sprintf("%0.2f",QscatR2);
        qalttext = sprintf("%0.0f",QscatR1);
        table = [table; [QS(j) QI(j) R1(j) R2(j) QscatR2 QscatR1]];
        if (i <= 0)
          text(R1(j)-0.05,R2(j),qtext,
               "horizontalalignment", "right",
               "verticalalignment", "middle",
               "fontname", "Times", "fontsize", 4);
          text(R1(j)+0.05,R2(j),"",#qalttext,
               "horizontalalignment", "left",
               "verticalalignment", "middle",
               "fontname", "Times", "fontsize", 4);
        end
      end
      #printf("Qsato Qint R1 R2 QscatR2 QscatR1\n");
      #table
      MeanQscatR2 = mean(QscatR2List);
      PctDevR2 = 100*std(QscatR2List)/MeanQscatR2;
      MeanQscatR1 = mean(QscatR1List);
      PctDevR1 = 100*std(QscatR1List)/MeanQscatR1;
      MeanQscatCombo = mean([QscatR2List, QscatR1List]);
      PctDevCombo = 100*std([QscatR2List, QscatR1List])/MeanQscatCombo;
      fmt = sprintf("%s %s %s %s\n",
                    "%s %8g;",
                    "%s = (%7g, %5.1f%%, %5.2f);",
                    "%s (%7g, %5.1f%%, %5.2f);",
                    "%s (%7g, %5.1f%%, %5.2f);");
      printf(fmt,
             "Sato Qscat:", uniqueQS(i), "Fehler Qscat: (avg, stdev, f/s)",
             MeanQscatR2, PctDevR2, MeanQscatR2/uniqueQS(i), "..R1..",
             MeanQscatR1, PctDevR1, MeanQscatR1/uniqueQS(i), "..avg..",
             MeanQscatCombo, PctDevCombo, MeanQscatCombo/uniqueQS(i)
            );
      ltext = sprintf("Q_{Fehler} = %5.0f +/-%3.0f%%", MeanQscatR2, PctDevR2);
      if (FehlerSatoRatio!=0)
        ltext = sprintf("Q_{Scat} = %5.0f",uniqueQS(i)*FehlerSatoRatio);
      end
      legend([h],ltext, "location", ScatLegLoc);
      legend("boxoff");
    end
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
    ltext = sprintf("Q_{Int} = %4.0f", uniqueQI(i));
    hl = legend([h],ltext, "location", IntLegLoc);
    set(hl, "fontsize", 7, "fontname", "Mono");
    legend("boxoff");
  end

  # Data Labels
  if (boolDataLabels)
    for (i=1:length(R1))
      text(R1(i), R2(i), sprintf("%d",i),
           "horizontalalignment", "right",
           "verticalalignment", "bottom",
           "fontname", "Courier", "fontsize", 3);
      text(R1(i), R2(i), sprintf("",B0(i)),
           "horizontalalignment", "left",
           "verticalalignment", "top",
           "fontname", "Times", "fontsize", 4);
    end
  end

  # Selected Q_scat computations:
  if (freq != 0)
    Rmult = 2.302585 * veloc / (2*pi*freq*(RangeFar-RangeNear));
  end

  # Labels and Titles
  # (Replace in fig scripts with: set(get(gca(),"title"),"string","New Title");
  #                               set(get(gca(),"xlabel"),"string","New Label");
  #  ...)
  labelfontname = "Helvetica";
  labelfontsize = 7;

  LabelX = sprintf("R1 = log10[E_{early} @ r = %g km / E_{late} @ r = %g km]",
                   RangeFar, RangeFar);
  LabelY = sprintf("R2 = log10[r^qE_{early} @ r = %g km / r^qE_{early} @ r = %g km]; q = %0.3f", RangeNear, RangeFar, GeoSpread);
  xlabel(LabelX,
        "fontweight", "bold", "fontsize", labelfontsize, "fontname", labelfontname);
  ylabel(LabelY,
        "fontweight", "bold", "fontsize", labelfontsize, "fontname", labelfontname);
  title("R_2 vs. R_1 in XXXXXX Model", "fontweight", "bold",
                    "fontsize", labelfontsize+3, "fontname", labelfontname);

end
