# normcurve_fitpowerlaw.m
#

# Given a NORMCURVE structure, find the powerlaw fit parameters (c and
# q) that best fit the curve over a selected range of elements (allows
# truncation of acharacteristic region).  The powerlaw fit function
# is:
#
#   Y = c * X .^ q
#
# Where Y is an energy series and X is a range (distance) series.
#
# The NORMCURVE structure contains two energy series, and soem
# metadata elements.  The elements of the energy series are
# characteristic energies of seismometers in a seismometer array. One
# energy series represents the integrated energies of the
# seismometers, and the other represents peak energies.  We are
# primarily interested in the integrated energies, but the peak
# energies may be of interest too, so we treat them both.  By
# smoothing the energy series via powerlaw fitting, we end up with a
# normalization curve that we can use to facilitate visual
# interpretation of traveltime curve amplitude in our image-density
# plots (as produced by arrayimage.m) and we have to ability to use a
# common normalization curve between multiple runs.
#
# NORMCURVE structure elements:
#
#  TimeWindow
#  RangeWindow
#  SummedEnergy
#  PeakEnergy
#
#  Added by this function:
#
#  PLCQ_Summed    // Power-Law C and Q, Summed trace
#  PLCQ_Peak      //  ''                Peak trace
#
function NC = normcurve_fitpowerlaw(NC, ibegin=1, iend=-1)

  # Defaults:
  if (iend==-1); iend = length(NC.SummedEnergy); end

  # Domain:
  X = linspace(NC.RangeWindow(1), NC.RangeWindow(2), length(NC.SummedEnergy)).';

  # Integrated Energy Curve:
  Y = log(NC.SummedEnergy(ibegin:iend));
  P = polyfit(log(X(ibegin:iend)), Y, 1);
  CQ = [exp(P(2)), P(1)];
  NC.SummedEnergy = CQ(1) * X .^ CQ(2);
  NC.PLCQ_Summed = CQ;

  # Peak Energy Curve:
  Y = log(NC.PeakEnergy(ibegin:iend));
  P = polyfit(log(X(ibegin:iend)), Y, 1);
  CQ = [exp(P(2)), P(1)];
  NC.PeakEnergy = CQ(1) * X .^ CQ(2);
  NC.PLCQ_Peak = CQ;

end
