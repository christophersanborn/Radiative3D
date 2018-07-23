# decimate.m
#
# Decimate a data series by picking max in min in a sliding window
# over a trace. This method of decimating is intended to preserve the
# visual envelope of a waveform when plotted, but is NOT intended to
# preserve any important analytical properties of the trace, like
# frequency makeup, etc.
#
# Given a window size W, decimate() will produce an output series
# approximately 2/W of the original size.  So for a 50% reduction,
# choose W=4.
#
function S = decimate(T,        # Input series, [t1 t2 t3 ... tn]
                      W         # Window width
                     )

  assert(isscalar(W));
  assert(W=floor(W));
  assert(W>=2);
  assert(isnumeric(T));
  assert(isvector(T));
  assert(length(T)>=W);

  N = floor(length(T)/W);
  R = mod(length(T),W);

  S=T(1:2);                     # Makes us a row or column vec to match T
  S(1:2) = maxmin(T(1:W));      # Get first pair
  for i = 1:(N-1)
    j = 2*i+1;
    k = W*i+1;
    S(j:j+1) = maxmin(T(k:k+W-1));
  end

  if (W>2 || R==1)        # What to do with the last bit needs to 
    S(end+1)= T(end);     # be thought out more. This will serve
  end                     # present purposes tho.

end

function mm = maxmin(V)
  [max, imax] = max(V);
  [min, imin] = min(V);
  ir = sort([imax,imin]);
  mm = V(ir);
end
