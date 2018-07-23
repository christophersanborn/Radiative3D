function annotate_array(velocs = [6.4,
                                  3.6,
                                  8.0,
                                  4.46],
                        phase = {"Pg",
                                 "Lg",
                                 "Pn",
                                 "Sn"},
                        range = 700,
                        linewidth = 1.33,
                        fontsize = 7.0)

if (min(size(velocs))==1)   # If just a list of velocities:
                            #
    L = sqrt(velocs.*velocs + 1);
    X = range * velocs./L;
    Y = range * (1.0)./L;

    for (j=1:length(velocs))
      x = X(j);
      y = Y(j);
      line([0 x], [0 y], "color", [0.3,0.3,0.3], "linewidth", linewidth);
      text(x, y, phase{j}, "fontsize",fontsize,
                           "color", [0.2,0.2,0.2],
                           "horizontalalignment", "left",
                           "verticalalignment", "middle");
    end


else    # Else assume velocs is an array with velocs, 
        #
    V  = velocs(:,1);  # Phase velocities
    T0 = velocs(:,2);  # Time axis intercepts
    R0 = velocs(:,3);  # Range on
    R1 = velocs(:,4);  # Range off

    M  = 1./V;         # Slope
    Y0 = T0 + M.*R0;
    Y1 = T0 + M.*R1;

    for (j=1:length(Y0))
      line([R0(j) R1(j)], [Y0(j) Y1(j)], 
                           "color", [0.3,0.3,0.3], "linewidth", linewidth);
      text(R1(j), Y1(j), phase{j}, 
                           "fontsize",fontsize,
                           "color", [0.2,0.2,0.2],
                           "horizontalalignment", "left",
                           "verticalalignment", "bottom");
    end


end ## endif
#####

end
#
