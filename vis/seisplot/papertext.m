# papertext.m
#
# Like the text() graphics function, but plots in paperspace rather
# than within an axes. Usage:
#
#   h = papertext(X, Y, "text", [prop, val], ... )
#
# X and Y both range from 0 to 1 and span the entire paperspace of the
# figure.
#
# Works by creating an invisible axes object that spans the whole
# figure, especially for paperspace annotations.  If such an axes
# object already exists, it will be reused instead of creating a new
# one.
#
function h = papertext(varargin)

  if (length(varargin) < 3)
     error("Usage: h = papertext(X, Y, \"text\", [prop, val], ... )");
  end

  hAx = gca();                    # Remember so we can restore later
  kids = get(gcf(), "children");  # Remember original child list too

  # Search for a paperspace axes, create if doesn't exist:
  hPaper = findobj(gcf(),"type", "axes", "position", [0 0 1 1], "visible", "off");
  if (length(hPaper)>0)
    hPaper = hPaper(end);
  else
    hPaper = axes('position',[0 0 1 1],'visible','off');
    kids = [kids; hPaper];        # Put new axis at bottom of z-order not top
    set(gcf(),"children",kids);   # (Fixes octave bug where invisible is
  end                             #  sometimes not honored.)

  # Draw text on paperspace
  axes(hPaper);
  h = text(varargin{});
  axes(hAx);

end
