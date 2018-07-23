## azimuth_deg.m
##
## Returns the azimuth in degrees east of north between two locations,
## nominally the location of a seismometer and the location of an
## event.
##
## Ultimately, this function will take into account the details of the
## Earth Coordinate system used to report locations, but for now we
## just assume the earth is flat, and that elevations don't matter.
## Thus we look just at the cartesian coordinates X and Y (ignoring
## Z), treating "north" as +y and "east" as +x.
##
## Inputs are assumed to be 1x3 matrices encoding [x, y, z].
##
##
function azimuth = azimuth_deg(dest, source)

  delta = dest - source;        # will get azi FROM source TO dest

  azimuth = atan2(delta(1), delta(2))*180/pi;
  if azimuth < 0
    azimuth += 360;     # favor 360* compass wheel
  end

end
