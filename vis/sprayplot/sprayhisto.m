# sprayhisto.m
#
function sprayhisto(ThetaPhiENUP,
                    ThetaPhiENUS,
                    Title = "Radiation by Azimuth",
                    NBins = 64
                   )

  AziP = (180/pi)*((pi/2)-ThetaPhiENUP(:,2));
  AziS = (180/pi)*((pi/2)-ThetaPhiENUS(:,2));

  clf();
  height = 3.50; # inches
  set(gcf(),"paperposition", [2.25 (11-height)/2 4.0 height]);
  set(gca(), "position", [0.11 0.15 0.84 0.74]);
  hold on;

  hist(AziS, NBins, "facecolor", "b");
  hist(AziP, NBins, "facecolor", "r");

  set(gca(), "xlim", [-90 270]);
  set(gca(), "xtick", linspace(-90,270,9));
  grid on;
  line([0 0], [0 axis()(4)], "linewidth", 3.0);

  title(Title);
  xlabel("Azimuth (degrees east of north)");

end
