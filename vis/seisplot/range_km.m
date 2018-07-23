## range_km.m
##
## Returns the "range" (distance along Earth's surface) between two
## locations, nominally the location of a seismometer and the location
## of an event.
##
## Ultimately, this function will take into account the details of the
## Earth Coordinate system used to report locations, but for now we
## just assume the earth is flat, and that elevations don't matter.
## Thus we just take the cartesian distance in X and Y (ignoring Z).
##
## Inputs dest and source are assumed to be 1x3 matrices encoding 
## [x, y, z].
##
##
function range = range_km(dest, source)

  delta = dest(:,1:2) - source(:,1:2);  # Get the x,y components only
  range = sqrt(delta * delta');

end
