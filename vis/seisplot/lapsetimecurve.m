## arrayimage.m
##
## Plot a lapse-time curve on existing plot axes
##
function Stats = lapsetimecurve (         # Stats = [B0, QScatS, QIntS, R1, R2]
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

  Stats = [Albedo, QScat, QInt, FehlerR1, FehlerR2];

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
  ## Deleted
end

function expl = getexplfrac(MT)

  tr = trace(MT);

  explpart = (tr^2) / 3;   ## Mag^2 of explosive part
  total = sum(sum(MT.*MT)); ## Mag^2 of whole tensor

  expl = explpart/total;
  if tr<0; expl*=-1; end

end
