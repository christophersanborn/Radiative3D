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
#
function seisplot (tracefile, gpow=1.0)

  SEIS = load(tracefile);
  if (isfield(SEIS,"SEIS")) SEIS=SEIS.SEIS; end # Catch old format
  if (isfield(SEIS,"AxesDesc"))
    AxesDesc = SEIS.AxesDesc;
  else
    AxesDesc = "XYZ";
  end

  # Look for MParams file:
  if (exist("out_mparams.octv","file")==2)
    MPAR=load("out_mparams.octv");
    NPhonCast = MPAR.NumPhonons;
  else
    NPhonCast = -1;
  end

  # Get unscaled signal maximum: (For plot range)
  SigMax_XYZ = max(max(SEIS.TraceXYZ));
  SigMax_PS  = max(max(SEIS.TracePS));
  SigMax     = max(SigMax_XYZ, SigMax_PS);

  # Get Energy and Count maxima:
  TotalE = sum(SEIS.TracePS,2);
  TotalC = sum(SEIS.CountPS,2);
  [EnergyMax,iEnergyMax] = max(TotalE);
  [CountMax,iCountMax] = max(TotalC);

  # Get Integrated Energy in each trace:
  dt = (SEIS.TimeWindow(2)-SEIS.TimeWindow(1)) / SEIS.NumBins;
  EIntXYZ = sum(SEIS.TraceXYZ,1) * dt;
  EIntPS  = sum(SEIS.TracePS,1)  * dt;
  CSumPS  = sum(SEIS.CountPS,1);

  # Apply gamma scaling:
  SEIS.TraceXYZ = SEIS.TraceXYZ .^ (gpow/2);
  SEIS.TracePS = SEIS.TracePS .^ (gpow/2);
  NormDetected = SigMax ^ (gpow/2);
  gpow_desc = sprintf("Scale: (amplitude)**%g", gpow);
  if gpow==1.0; gpow_desc="Scale: Amplitude"; end
  if gpow==2.0; gpow_desc="Scale: Energy"; end

  # Apply Normalization:
  NormGoal = 190;
  Staging  = 200;
  NormFactor = NormGoal / NormDetected;
  NormFactPS = NormGoal / max(max(SEIS.CountPS));

  EEX = NormFactor * SEIS.TraceXYZ(:,1);  ## Envelope Energy X
  EEY = NormFactor * SEIS.TraceXYZ(:,2);  ## Envelope Energy Y
  EEZ = NormFactor * SEIS.TraceXYZ(:,3);  ## Envelope Energy Z
  EEP = NormFactor * SEIS.TracePS(:,1);   ## Envelope Energy P
  EES = NormFactor * SEIS.TracePS(:,2);   ## Envelope Energy S
  ECP = NormFactPS * SEIS.CountPS(:,1);   ## Envelope Phonon-Count P
  ECS = NormFactPS * SEIS.CountPS(:,2);   ## Envelope Phonon-Count S


  ##
  ## Plot Prep: Layout and Style
  ##

  linethin   = 1.5;
  linethick1 = 3.75;
  fonttitle  = 9.5;
  fontxylabel = 9;
  fontaxes  = 9;
  fontinset = 9;
  fonttiny  = 7;
  color_maxmark = [1.0,0.6,0.2];

  figinit(5.0, 3.75,                # Clear and initialize 
       "paperunits", "inches",      # a 5.0" x 3.75" figure
       "visible", "off");
  axes("fontweight", "bold",        # Create and setup an axes object
       "fontsize", fontaxes,        #
       "linewidth", linethin, "ytick", [5 4 3 2 1]*Staging,
       "yticklabel", {AxesDesc(1), AxesDesc(2), AxesDesc(3), "P", "S"},
       "nextplot", "add");  # equiv to "hold on"

  Flatline = zeros(length(EEX),1);
  x_min = SEIS.TimeWindow(1);
  x_max = SEIS.TimeWindow(2);
  x_data = linspace(x_min, x_max, length(EEX));
  axis([x_min x_max 0.05*Staging, 6.35*Staging]);
  xlabel("Time (s)", "fontsize", fontxylabel);

  ## 
  ## Begin Plotting:
  ##

  linewidth = linethick1;

  ## P-component Count Plot:
  Baseline = Flatline + 2*Staging;
  color = [1.00 0.70 0.90];     # count color
  fillplot(x_data',Baseline+ECP,color,2*Staging);

  ## S-component Count Plot:
  Baseline = Flatline + 1*Staging;
  color = [0.70 0.90 1.00];     # count color
  fillplot(x_data',Baseline+ECS,color,1*Staging);

  ## Vertical Max Lines
  line(x_data(iEnergyMax)*[1,1], Staging*[0.8,5.7], 
       "linestyle", ":", "linewidth", linethin, "color", color_maxmark);
  line(x_data(iCountMax)*[1,1], Staging*[0.8,5.7], 
       "linestyle", ":", "linewidth", linethin, "color", color_maxmark);

  ## X-component Plot:
  Baseline = Flatline + 5*Staging;
  color = [0.85 0.00 0.00];
  plot(x_data, Baseline, "linewidth", linethin);
  set(get(gca(),"children")(1),"color", color);
  plot(x_data, Baseline+EEX, 'LineWidth', linewidth);
  set(get(gca(),"children")(1),"color", color);

  ## Y-component Plot:
  Baseline = Flatline + 4*Staging;
  color = [0.00 0.60 0.00];
  plot(x_data, Baseline, "linewidth", linethin);
  set(get(gca(),"children")(1),"color", color);
  plot(x_data, Baseline+EEY, 'LineWidth', linewidth);
  set(get(gca(),"children")(1),"color", color);

  ## Z-component Plot:
  Baseline = Flatline + 3*Staging;
  color = [0.00 0.10 0.75];
  plot(x_data, Baseline, "linewidth", linethin);
  set(get(gca(),"children")(1),"color", color);
  plot(x_data, Baseline+EEZ, 'LineWidth', linewidth);
  set(get(gca(),"children")(1),"color", color);

  ## P-component Plot:
  Baseline = Flatline + 2*Staging;
  color = [0.90 0.00 0.60];     # energy/amplitude color
  plot(x_data, Baseline, "linewidth", linethin);
  set(get(gca(),"children")(1),"color", color);
  plot(x_data, Baseline+EEP, 'LineWidth', linewidth);
  set(get(gca(),"children")(1),"color", color);

  ## S-component Plot:
  Baseline = Flatline + 1*Staging;
  color = [0.00 0.70 0.90];     # energy/amplitude color
  plot(x_data, Baseline, "linewidth", linethin);
  set(get(gca(),"children")(1),"color", color);
  plot(x_data, Baseline+EES, 'LineWidth', linewidth);
  set(get(gca(),"children")(1),"color", color);

  ## S count plot:
  
  ##
  ## Text Annotations:
  ##

  # (First convert data from Cartesian to azimuth)
  RSeisLoc = SEIS.Location - SEIS.EventLoc; # Location relative to event
  azimuth = atan2(RSeisLoc(1), RSeisLoc(2))*180/pi;
  if azimuth < 0
    azimuth += 360;
  end
  Range = sqrt(RSeisLoc(1)^2 + RSeisLoc(2)^2);

  # Precompute width of a time bin in cycles:
  tbinwidth = (SEIS.Frequency * dt);

  # Titles and axes labels:
  title(sprintf("Location: %0.2f km at %0.2f deg az from source", 
                Range, azimuth),
        "fontsize", fonttitle);

  fontfam = "Courier";

  # Curve integrals:
  fontspec = {"fontsize", fonttiny, "horizontalalignment", "right",...
              "verticalalignment", "bottom", "fontname", fontfam};
  htxt = {};
  htxt{end+1} = text(x_max*0.99, 2*Staging+110,
                     sprintf("Sum(E*dt) X: %8.2e", EIntXYZ(1))); 
  htxt{end+1} = text(x_max*0.99, 2*Staging+70,
                     sprintf("Y: %8.2e", EIntXYZ(2)));
  htxt{end+1} = text(x_max*0.99, 2*Staging+30,
                     sprintf("Z: %8.2e", EIntXYZ(3)));
  htxt{end+1} = text(x_max*0.99, 1*Staging+110,
                     sprintf("Sum(E*dt) P: %8.2e", EIntPS(1))); 
  htxt{end+1} = text(x_max*0.99, 1*Staging+70,
                     sprintf("S: %8.2e", EIntPS(2)));
  htxt{end+1} = text(x_max*0.99, 1*Staging+30,
                     sprintf("P+S: %8.2e", sum(EIntPS)));
  set([htxt{}],fontspec{});

  # Count statistics:
  text_phoncaught = sprintf("Total Caught: %8.2e", sum(CSumPS));
  if (NPhonCast > 0) 
    text_phoncast = sprintf("Phonons Cast: %8.2e", NPhonCast);
    text_catchrate = sprintf("Catch Rate: %8.2e", sum(CSumPS)/NPhonCast);
  else 
    text_phoncast = "Phonons Cast: [unknown]";
    text_catchrate = "";
  end
  fontspec = {"fontsize", fonttiny, "horizontalalignment", "right",...
              "verticalalignment", "top", "fontname", fontfam};
  htxt = {};
  htxt{end+1} = text(x_max*0.99, 4*Staging+190, text_phoncast);
  htxt{end+1} = text(x_max*0.99, 4*Staging+150, text_phoncaught);
  htxt{end+1} = text(x_max*0.99, 4*Staging+110, text_catchrate);
  set([htxt{}],fontspec{});

  # Additional annotations:
  fontspec = {"fontname", fontfam};
  htxt = {};
  htxt{end+1} = text(SEIS.TimeWindow(2)*.05,6*Staging+20, 
       sprintf("Energy Max: %8.2e at %0.2fs", EnergyMax, x_data(iEnergyMax)),
       "fontsize", fontinset);
  htxt{end+1} = text(SEIS.TimeWindow(2)*.05,6*Staging-30, 
       sprintf("Count Max:  %-8d at %0.2fs", CountMax, x_data(iCountMax)),
       "fontsize", fontinset);
  htxt{end+1} = text(SEIS.TimeWindow(2)*.65,6*Staging+20, 
       gpow_desc, "horizontalalignment", "left",
       "fontsize", fontinset);
  htxt{end+1} = text(SEIS.TimeWindow(2)*.65,6*Staging-30,
       sprintf("Time Bin: %0.2f s", tbinwidth/SEIS.Frequency),
       "horizontalalignment", "left",
       "fontsize", fontinset);
  htxt{end+1} = text(SEIS.TimeWindow(2)*.05,140, 
       sprintf("Frequency = %0.2f Hz", SEIS.Frequency),
       "fontsize", fontinset);
  explfrac = getexplfrac(SEIS.EventMT);
  htxt{end+1} = text(SEIS.TimeWindow(2)*.05,95, 
       sprintf("Event Isotropy: %0.0f%%", explfrac*100),
       "fontsize", fontinset);
  htxt{end+1} = text(SEIS.TimeWindow(2)*.5,140, 
       sprintf("Gather Radius P = %0.1f km", SEIS.GatherRadius(1,2)),
       "fontsize", fontinset);
  htxt{end+1} = text(SEIS.TimeWindow(2)*.5,95, 
       sprintf("Gather Radius S = %0.1f km", SEIS.GatherRadius(2,2)),
       "fontsize", fontinset);
  set([htxt{}],fontspec{});
  grid("off");

end

function h = fillplot(X,Y,C,Y0=0)

  h = fill([X(1); X; X(end)],
           [  Y0; Y; Y0    ], C);
  set(h, "linestyle", "none", "edgecolor", C);

end

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
