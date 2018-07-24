## baseplot_gridWCGrangeelev.m
##
## Plots grid geometry range and elevation mesh onto an existing axis.
##
function baseplot_gridWCGrangeelev(
           GG,            # The grid, indexed in iRange,iAzi,iZ,iAttr
           ColorMap = 0,  # If Nx3 matrix, then we fill plot instead
                          #  of mesh plot, and each row of ColorMap
                          #  colors a "layer" in the model.
           RotMat = [1 0; 0 1]  # Optional rotation matrix, e.g. to adjust
                                # orientation of an OCS RAW grid.
         )

  [nx,ny,nz,nattr] = size(GG);
  fillplot = false;
  if (size(ColorMap,2)==3)
    fillplot = true;
  end

  # Make Range-Elevation arrays:
  R = zeros(nx,nz);
  E = zeros(nx,nz);
  for ix = 1:nx
    for iz = 1:nz
      x = GG(ix,1,iz,1);
      y = GG(ix,1,iz,2);
      z = GG(ix,1,iz,3);
      sg = sign((pi/2)-abs(atan2(y,x)));  # -1 if "behind" origin wrt due east
      R(ix,iz) = sqrt(x^2+y^2) * sg;      #
      E(ix,iz) = z;
      rotXY = RotMat*[R(ix,iz);E(ix,iz)]; # Transform (e.g. rotate) node points
      R(ix,iz) = rotXY(1);                # about origin point.
      E(ix,iz) = rotXY(2);                #
    end
  end

  # Plot polygons to represent WCG unit cells:
  for ir = 1:(nx-1)
    for ie = 1:(nz-1)
      X = [R(ir,ie), R(ir+1,ie), R(ir+1,ie+1), R(ir,ie+1), R(ir,ie)];
      Z = [E(ir,ie), E(ir+1,ie), E(ir+1,ie+1), E(ir,ie+1), E(ir,ie)];
      if (fillplot)
        Color = ColorMap(min(ie,size(ColorMap,1)),:);
        hpatch = fill(X,Z,Color);
        edgeColor = Color * 0.25;
        set(hpatch,"edgecolor",edgeColor);
        set(hpatch,"linewidth",0.5);
      else
        line (X,Z,"linewidth",1.25,"color",[0,0,0]);
      end
    end
  end

end
