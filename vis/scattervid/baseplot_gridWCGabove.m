## baseplot_gridWCGabove.m
##
## Plots grid geometry map-view mesh onto an existing axis.
##
function baseplot_gridWCGabove(
           GG,            # The grid, indexed in iRange,iAzi,iZ,iAttr
           ColorMap = 0   # If Nx3 matrix, then we fill plot instead
                          #  of mesh plot, and each row of ColorMap
                          #  colors a vertical strip in the model.
         )

  [nr,na,nz,nattr] = size(GG);
  fillplot = false;
  if (size(ColorMap,2)==3)
    fillplot = true;
  end

  # Make Easting, Northing arrays:
  E = zeros(nr,na);
  N = zeros(nr,na);
  for ir = 1:nr
    for ia = 1:na
      x = GG(ir,ia,1,1);
      y = GG(ir,ia,1,2);
      E(ir,ia) = x;
      N(ir,ia) = y;
    end
  end

  # Plot polygons to represent WCG unit cells:
  for ir = 1:(nr-1)
    for ia = 1:(na-1)
      X = [E(ir,ia), E(ir+1,ia), E(ir+1,ia+1), E(ir,ia+1), E(ir,ia)];
      Y = [N(ir,ia), N(ir+1,ia), N(ir+1,ia+1), N(ir,ia+1), N(ir,ia)];
      if (fillplot)
        Color = ColorMap(min(ir,size(ColorMap,1)),:);
        fill(X,Y,Color);
      else
        line (X,Y,"linewidth",1.5,"color",[0,0,0]);
      end
    end
  end

end
