## read_gridgeom.m
##
## Reads a ascii dump of the model grid compacts the geometrical parts
## (XYZ coords) into a matrix that can be used for grid plotting, etc.
##
## Matrix is called GG for grid geometry.  This may be complemented by
## an additional matrix GA containing the grid attributes (velocities,
## etc.), but at present this file does not deal with that.
##
function GG = read_gridgeom(
                gridfile        # name of file with grid data
              )

  GG_Raw = load(gridfile);

  ibase = 0;                            # (Assume zero, but might not be.)
  NX = (1-ibase)+max(GG_Raw,[],1)(1);   # Grid dimensions
  NY = (1-ibase)+max(GG_Raw,[],1)(2);   #
  NZ = (1-ibase)+max(GG_Raw,[],1)(3);   #

  NAttrDefs = size(GG_Raw,1);

  GG = zeros(NX, NY, NZ, 3);            # 12 Attributes (XYZ, VpVs,
                                        # Dens, QsQp, NEAK) in 3 indices,
                                        # but we are here only intinding 
                                        # to retain three of them (X,Y,Z)

  for (j=1:NAttrDefs)                   # Pack each line of raw grid into 
      ix = (1-ibase)+GG_Raw(j,1);       # appropriate cubby holes.  Repeats
      iy = (1-ibase)+GG_Raw(j,2);       # will overwrite, but this is OK, as
      iz = (1-ibase)+GG_Raw(j,3);       # the XYZ's shouldn't be different.
      for (k=(1:3))
        GG(ix,iy,iz,k) = GG_Raw(j,3+k);
      end
  end

end
