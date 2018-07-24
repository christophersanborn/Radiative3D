#
#
function lapsetimebaseplot()

  figinit(6.0, 6.0,                # Clear and initialize 
       "paperunits", "inches",      # a 5.0" x 3.75" figure
       "visible", "on");           #
  axes("fontsize", 12,        # Create and initialize axes obj
       "linewidth", 2, "nextplot", "add");

  global linestyle;
  linestyle.width = 3;
  linestyle.style = "-";

  axis([0 250, 0.01 100]);
  title("Blue: Early window (x to x s); Green: Late window (x to x s)",
        "fontweight", "bold");
  ylabel("r^2E_{1,2}(r) / r^2E_1(0)", "fontweight", "Bold");
  xlabel ("Distance (km)", "fontweight", "Bold");

end
