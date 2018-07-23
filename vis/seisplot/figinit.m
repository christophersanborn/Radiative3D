# figinit.m
#
# Initialize a figure in a consistent way, setting the size of the
# figure in inches, and other properties if desired.
#
# Usage: figinit(width, height, [property], [value], ...)
#
function figinit(varargin)

  [reg, PV] = parseparams(varargin);

  if (length(reg)!=2)   # Want exactly two "regular" args (width & height)
     error("Usage: [out] = figinit(width, height, property, value, ...)");
  end

  width = reg{1};       # Assumed inches, but can override by passing 
  height = reg{2};      #   "paperunits", "centimeters", ...

  graphics_toolkit "gnuplot";   # Otherwise it crashes when run over
                                # SSH connection. (Affects newer Octave
                                # versions with OpenGL and Qt support.)

  clf();                # Clear it all...

  if (length(PV)>1)     # Set all requested properties. (if units
    set(gcf(), PV{});   # changed, will get papersize in correct units 
  end                   # before we query it)

                        # BTW, not sure what happens if user sets 
                        # paperorientation to landscape...

  papersize = get(gcf(), "papersize");  # Get papersize to compute margins
  pwidth = papersize(1);
  pheight = papersize(2);
  pwmargin = (pwidth-width)/2;
  phmargin = (pheight-height)/2;

  paperposition = [pwmargin, phmargin, width, height];

  set(gcf(), "paperposition", paperposition);

endfunction
