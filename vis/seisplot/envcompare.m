# envcompare.m
#
# Plots two envelopes in for comparison.  One could be synthetic, one
# data, e.g.
#
# Plots using fill-to-axes plots, one in front of the other, both
# above and below axis, with below axis switching the order.
#
# 
function envcompare(A1,  # Series 1, typically synthetic amplitude env
                    A2,  # Series 2, typically real-data env
                    twindow, # Time window, e.g, [0 600]
                    rscale,  # A2 norm multiplier (top,bottom)
                    smooth=[8,0], # Num times to smooth date (A1, A2)
         titletext="Station MAK channel BHZ: envelope comparison at 3.0 Hz"
                   )

  x1=linspace(twindow(1), twindow(2), length(A1))';
  x2=linspace(twindow(1), twindow(2), length(A2))';

  # Normalize and Smooth:
  for i = 1:smooth(1)
    A1 = ThreePSmooth(A1);
  end
  for i = 1:smooth(2)
    A2 = ThreePSmooth(A2);
  end
  A1 ./= max(A1);
  A2 ./= max(A2);

  # Initialize plot
  clf();
  set(gcf(),"paperposition", [0.25 2.5 8.0 3.0]); # double-wide
  hold on;

  # Colors:
  A1fillc = [.3 .6 .8];
  A2fillc = [.8 .6 .2];
  A1linec = 0.7*(A1fillc.^2);
  A2linec = 0.7*(A2fillc.^2);
  shadowfact = 1.0;#0.5;

  # Rear:
  h1 = fillplot(x2, rscale(1)*A2, shadowfact*A2fillc);
  h2 = fillplot(x1,-A1, shadowfact*A1fillc);
  # Rear-ish
  plot (x2,rscale(1)*A2,"linewidth", 4.0, "color", A2linec);
  plot (x1,-A1,"linewidth", 4.0, "color", A1linec);
  # Front:
  h3 = fillplot(x1, A1, A1fillc);
  h4 = fillplot(x2,-rscale(2)*A2, A2fillc);
  # Front-ish
  hA1 = plot (x1,A1,"linewidth", 4.0, "color", A1linec);
  hA2 = plot (x2,-rscale(2)*A2,"linewidth", 4.0, "color", A2linec);
  set(h1, "facealpha", 0.5);
  set(h2, "facealpha", 0.5);
  set(h3, "facealpha", 0.5);
  set(h4, "facealpha", 0.5);

  xlabel("Time (s)", "fontweight", "bold");
  ylabel("Env Amplitude", "fontweight", "bold");
  title(titletext, "fontweight", "bold");
  set(gca(),"fontweight", "bold");
  axis([1 1 1.2 1.2].*axis(),"ticx");
  hL = legend([hA1, hA2], "Synthetic", "Data", "location", "northeast");
  set(hL,"fontweight","bold","linewidth",800.0); # linewidth doesnt work...
  set(get(hL,"children"),"linewidth",800.0);     # ''

end

function h = fillplot(X,Y,C)

  h = fill([X(1); X; X(end)],
           [   0; Y; 0     ], C);
  set(h, "linestyle","none", "edgecolor", C);

end

function S=ThreePSmooth(D)
  S = zeros(length(D),1);
  S(1) = .75*D(1) + 0.25*D(2);
  S(end) = .75*D(end) + 0.25*D(end-1);
  for i = 2:(length(D)-1)
      S(i) = 0.25*D(i-1) + 0.5*D(i) + 0.25*D(i+1);
  end
end
